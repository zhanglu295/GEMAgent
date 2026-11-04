#!/usr/bin/env python3
"""LLM adjudication for corrected pure-OR GPR candidates.

Deterministic EC grouping is supplied as evidence only. The LLM decides whether
each OR branch is biologically plausible for the reaction. No GE truth or
curated-model information enters the prompt.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
import sys
import time
import urllib.request
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


SYSTEM = """You are a senior metabolic model curator.

Your task is to review one pure-OR gene-protein-reaction rule. Every listed gene
is independently sufficient for the reaction. Decide which branches are
biologically plausible for the specific reaction and which are likely
mis-attached by automatic reconstruction.

Use the reaction chemistry, reaction name, seed EC, KOfam KO definitions, gene
EC assignments, transporter specificity, known enzyme multifunctionality, and
the logical structure of the current OR list. You may infer a mismatch from
converging annotation and biochemical evidence even when no source explicitly
says "this gene is mis-attached".

Decision policy:
- KEEP is the default when the branch is plausible or evidence is insufficient.
- REMOVE only when the annotation or biochemical specificity positively
  conflicts with this reaction and at least one other branch remains plausible.
- Never remove every gene in the OR rule.
- Do not edit the reaction itself. Do not use essentiality labels, curated GPRs,
  or any hidden benchmark answer. They are not supplied and must not be guessed.
- Distinguish a different EC in the same broad reaction family from an unrelated
  function. A shared first two EC levels can support a plausible alternative,
  while transport substrate class and reaction chemistry must still agree.

Return one decision per gene with a concise, auditable rationale.
"""


OUTPUT_SPEC = """Return ONLY valid JSON:
{
  "decisions": [
    {
      "gene": "<exact model gene id>",
      "verdict": "KEEP|REMOVE|UNCERTAIN",
      "confidence": 0.00,
      "reason": "<concise biochemical or annotation-based rationale>",
      "evidence_basis": "annotation_chemistry|transport_specificity|isozyme_logic|model_structure|uncertain"
    }
  ]
}
Include every listed gene exactly once. Do not output markdown.
"""


def norm_rid(value: str) -> str:
    value = str(value).strip()
    return value[2:] if value.startswith("R_") else value


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_json_response(text: str) -> list[dict]:
    cleaned = str(text or "").strip()
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned)
        cleaned = re.sub(r"\s*```$", "", cleaned)
    data = json.loads(cleaned)
    decisions = data.get("decisions")
    if not isinstance(decisions, list):
        raise ValueError("response has no decisions list")
    return decisions


def load_kofam_definitions(path: Path) -> dict[str, list[str]]:
    out = defaultdict(list)
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 6)
            if len(parts) < 7 or parts[0] != "*":
                continue
            gene = parts[1].strip()
            ko = parts[2].strip()
            definition = parts[6].strip().strip('"')
            out[gene].append(f"{ko}: {definition}")
    return out


def build_prompt(row: dict, genes: list[dict], manifest: dict[str, dict], organism: str) -> str:
    rid = row["reaction_id"]
    context = manifest.get(norm_rid(rid), {})
    name = context.get("name") or row.get("name") or "(reaction name unavailable)"
    seed_ec = context.get("ec") or "(not supplied)"
    clusters = row.get("clusters") or "(no full-EC consensus)"

    lines = [
        SYSTEM,
        "",
        "## Reaction context",
        f"- Organism: {organism}",
        f"- Reaction ID: {rid}",
        f"- Reaction name: {name}",
        f"- SEED EC: {seed_ec}",
        f"- Current OR branch count: {row['n_genes']}",
        f"- EC annotation groups in the current OR list: {clusters}",
        "",
        "## Current gene branches",
    ]
    for idx, gene in enumerate(genes, 1):
        kos = gene.get("kos") or "no KOfam hit"
        ecs = gene.get("ecs") or "no EC"
        label = gene.get("gene_label") or ""
        label_text = f" label={label}" if label else ""
        definitions = "; ".join(gene.get("ko_definitions", [])) or "no KOfam definition"
        lines.append(
            f"{idx}. model_gene={gene['gene']}{label_text}; "
            f"KOfam_gene={gene.get('kofam_gene') or 'unmapped'}; "
            f"KO={kos}; gene_EC={ecs}; definitions={definitions}"
        )
    lines.extend([
        "",
        "## Task",
        "For each branch, judge whether its encoded function can plausibly catalyze "
        "this exact reaction or transport process. Use the fact that this is a pure "
        "OR rule: two legitimate isozymes may share the same EC, while a branch with "
        "a different enzyme family or incompatible substrate class is likely "
        "mis-attached. If the reaction name is generic or annotations conflict, use "
        "the dominant annotation group only as supporting evidence, not as truth.",
        "",
        OUTPUT_SPEC,
    ])
    return "\n".join(lines)


def validate_decisions(raw: list[dict], expected_genes: list[str]) -> list[dict]:
    by_gene = {}
    for item in raw:
        if not isinstance(item, dict):
            continue
        gene = str(item.get("gene", "")).strip()
        if gene not in expected_genes:
            continue
        verdict = str(item.get("verdict", "")).strip().upper()
        if verdict not in {"KEEP", "REMOVE", "UNCERTAIN"}:
            verdict = "UNCERTAIN"
        try:
            confidence = float(item.get("confidence", 0.0))
        except (TypeError, ValueError):
            confidence = 0.0
        by_gene[gene] = {
            "gene": gene,
            "verdict": verdict,
            "confidence": max(0.0, min(1.0, confidence)),
            "reason": str(item.get("reason", "") or "").strip(),
            "evidence_basis": str(item.get("evidence_basis", "") or "uncertain").strip(),
        }

    out = []
    for gene in expected_genes:
        out.append(by_gene.get(gene, {
            "gene": gene,
            "verdict": "UNCERTAIN",
            "confidence": 0.0,
            "reason": "LLM omitted this branch",
            "evidence_basis": "uncertain",
        }))

    if all(d["verdict"] == "REMOVE" for d in out):
        for d in out:
            d["verdict"] = "UNCERTAIN"
            d["reason"] = "Safety fallback: an OR rule cannot be fully removed."
    return out


_LLM_ENV = None


def llm_env():
    global _LLM_ENV
    if _LLM_ENV is not None:
        return _LLM_ENV
    env = {}
    env_path = os.getenv("LLM_ENV_FILE") or str(Path(__file__).resolve().parents[5] / ".env")
    with open(env_path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, value = line.split("=", 1)
                env[key] = value.strip().strip('"').strip("'")
    _LLM_ENV = env
    return env


def call_llm(messages):
    env = llm_env()
    url = env["LLM_BASE_URL"].rstrip("/") + "/chat/completions"
    payload = {
        "model": env.get("LLM_MODEL_THINK") or env.get("LLM_MODEL", "deepseek-v4-flash"),
        "messages": messages,
        "temperature": 0,
        "max_tokens": 4096,
        "response_format": {"type": "json_object"},
        "thinking": {"type": "enabled"},
    }
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Authorization": "Bearer " + env["LLM_API_KEY"],
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(request, timeout=120) as response:
        data = json.loads(response.read().decode("utf-8"))
    return data["choices"][0]["message"]["content"] or ""


def review_one(row, gene_rows, manifest, prompt_dir, max_retries, organism):
    rid = row["reaction_id"]
    genes = gene_rows[rid]
    prompt = build_prompt(row, genes, manifest, organism)
    prompt_path = prompt_dir / f"{rid.replace('/', '_')}.txt"
    prompt_path.write_text(prompt, encoding="utf-8")
    messages = [{"role": "user", "content": prompt}]
    last_error = ""
    raw_text = ""
    review_failed = True
    for attempt in range(max_retries + 1):
        try:
            raw_text = call_llm(messages)
            decisions = validate_decisions(
                parse_json_response(raw_text),
                [g["gene"] for g in genes],
            )
            return rid, decisions, raw_text, False
        except Exception as exc:
            last_error = str(exc)
            if attempt < max_retries:
                time.sleep(2 ** attempt)
    decisions = [
        {
            "gene": g["gene"],
            "verdict": "UNCERTAIN",
            "confidence": 0.0,
            "reason": f"LLM error: {last_error[:240]}",
            "evidence_basis": "uncertain",
        }
        for g in genes
    ]
    return rid, decisions, raw_text, review_failed


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit_dir", required=True)
    ap.add_argument("--legacy_manifest", required=True)
    ap.add_argument("--kofam", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--organism", required=True)
    ap.add_argument("--cohort", choices=("pilot", "review", "all"), default="pilot")
    ap.add_argument("--limit", type=int, default=30)
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--max_retries", type=int, default=2)
    ap.add_argument("--llm_config_dir", required=True)
    args = ap.parse_args()

    if args.workers < 1:
        ap.error("--workers must be >= 1")
    if args.limit < 0:
        ap.error("--limit must be >= 0")

    audit_dir = Path(args.audit_dir)
    out_dir = Path(args.out_dir)
    prompt_dir = out_dir / "prompts"
    out_dir.mkdir(parents=True, exist_ok=True)
    prompt_dir.mkdir(parents=True, exist_ok=True)
    sys.path.insert(0, args.llm_config_dir)

    reactions = read_tsv(audit_dir / "or2_reactions.tsv")
    gene_rows_flat = read_tsv(audit_dir / "or2_genes.tsv")
    kofam_definitions = load_kofam_definitions(Path(args.kofam))
    manifest = {
        row["reaction_id"]: row
        for row in read_tsv(Path(args.legacy_manifest))
    }
    gene_rows = defaultdict(list)
    for gene in gene_rows_flat:
        gene["ko_definitions"] = kofam_definitions.get(gene.get("kofam_gene", ""), [])
        gene_rows[gene["reaction_id"]].append(gene)

    if args.cohort == "review":
        cohort = [r for r in reactions if r["category"].endswith("needs_review")]
    elif args.cohort == "all":
        cohort = list(reactions)
    else:
        candidates = [
            r for r in reactions
            if r["category"] == "confident_removal_candidate"
        ]
        review = [
            r for r in reactions
            if r["category"].endswith("needs_review")
        ]
        review.sort(
            key=lambda r: (
                -int(r["n_genes"]),
                -int(r["n_ec_annotated"]),
                r["reaction_id"],
            )
        )
        cohort = candidates + review

    seen = set()
    selected = []
    for row in cohort:
        if row["reaction_id"] in seen:
            continue
        seen.add(row["reaction_id"])
        selected.append(row)
        if args.limit and len(selected) >= args.limit:
            break

    if selected:
        try:
            llm_env()
        except Exception as exc:
            print(f"OR2 LLM configuration error: {exc}", file=sys.stderr)
            return 2

    results = {}
    raw_responses = {}
    failed_reviews = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(
                review_one,
                row,
                gene_rows,
                manifest,
                prompt_dir,
                args.max_retries,
                args.organism,
            ): row
            for row in selected
        }
        completed = 0
        for future in as_completed(futures):
            row = futures[future]
            rid, decisions, raw_text, failed_review = future.result()
            if failed_review:
                failed_reviews += 1
            results[rid] = decisions
            raw_responses[rid] = raw_text
            completed += 1
            if completed % 5 == 0 or completed == len(selected):
                print(f"progress {completed}/{len(selected)}", flush=True)

    decisions_path = out_dir / "llm_decisions.tsv"
    with decisions_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "reaction_id", "gene", "verdict", "confidence", "evidence_basis", "reason"
        ])
        for row in selected:
            rid = row["reaction_id"]
            for decision in results.get(rid, []):
                writer.writerow([
                    rid,
                    decision["gene"],
                    decision["verdict"],
                    f"{decision['confidence']:.4f}",
                    decision["evidence_basis"],
                    decision["reason"],
                ])

    (out_dir / "raw_responses.json").write_text(
        json.dumps(raw_responses, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    llm_model = ""
    if selected:
        llm_model = (
            llm_env().get("LLM_MODEL_THINK")
            or llm_env().get("LLM_MODEL", "")
        )

    run_manifest = {
        "cohort": args.cohort,
        "limit": args.limit,
        "workers": args.workers,
        "llm_model": llm_model,
        "organism": args.organism,
        "reviewer_code_sha256": sha256(Path(__file__).resolve()),
        "source_reactions": str(audit_dir / "or2_reactions.tsv"),
        "source_reactions_sha256": sha256(audit_dir / "or2_reactions.tsv"),
        "source_genes": str(audit_dir / "or2_genes.tsv"),
        "source_genes_sha256": sha256(audit_dir / "or2_genes.tsv"),
        "n_selected_reactions": len(selected),
        "failed_reviews": failed_reviews,
        "selected_reactions": [r["reaction_id"] for r in selected],
        "truth_firewall": "pass: only model/KOfam annotations and reaction context supplied",
        "decision_file_sha256": sha256(decisions_path),
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    counts = defaultdict(int)
    for decisions in results.values():
        for decision in decisions:
            counts[decision["verdict"]] += 1
    print(json.dumps({
        "reactions": len(selected),
        "failed_reviews": failed_reviews,
        "verdict_counts": dict(sorted(counts.items())),
        "out_dir": str(out_dir),
    }, ensure_ascii=False, indent=2))
    if selected and failed_reviews == len(selected):
        print("OR2 review failed for every selected reaction.", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
