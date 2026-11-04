#!/usr/bin/env python3
"""
GEMAgent — post-gapfill LLM bypass review (iterative).

Reads {genome}-gapfill-added.json (reactions that gapfilling added), asks the
LLM to judge whether each is a legitimate reaction for the target organism or a
wrong pathway / bypass leak that would mask gene essentiality.  Reactions the
LLM excludes are appended to a per-organism exclusion list and gapfilling is
re-run; the loop continues until the LLM stops excluding reactions or the max
iteration cap is reached.

This never re-runs the candidate-review LLM. It is a separate, gapfill-specific
review whose INCLUDE/EXCLUDE semantics target bypass detection, not pathway
evidence accumulation.

OPT-2 (confidence layering): only HIGH-confidence EXCLUDE decisions are appended
to the permanent exclusion list. MEDIUM/LOW-confidence excludes are recorded in
the decisions JSON for transparency but do NOT permanently remove the reaction,
so a geneless-but-legitimate reaction the LLM is unsure about stays in the model.

OPT-3 (pathway context): the prompt is enriched with the organism's predicted
pathway inventory ({genome}-all-Pathways.tbl), so GENELESS reactions (no BLAST
hit) are judged by PATHWAY MEMBERSHIP rather than by the absence of gene
evidence. Absence of a BLAST hit is usually an annotation gap, not a foreign
pathway.
"""

from __future__ import annotations

import argparse
import glob
import html
import json
import os
import re
import subprocess
import sys
import time


current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)
sys.path.insert(0, os.path.dirname(current_dir))   # <project>/code, for paths

from llm_config import llm
from seed_db import _load_seed_db
import paths
import taxonomy

RSCRIPT = paths.RSCRIPT
CODE_DIR = paths.CODE_DIR

GAPFILL_SYSTEM_PROMPT = """You are a senior metabolic network reconstruction specialist.
You are reviewing reactions that a GAPFILLING algorithm added to a genome-scale
metabolic model. Gapfilling adds reactions ONLY because FBA requires them to
produce biomass — it has NO evidence that the organism actually possesses them.

Your job: decide for each added reaction whether it is a REAL metabolic step of
the target organism (INCLUDE) or a WRONG-PATHWAY / BYPASS LEAK (EXCLUDE).

## EXCLUDE when any of these apply

1. WRONG PATHWAY: the reaction belongs to a pathway the organism is KNOWN NOT to
   use. Example: a C4/Shemin 5-aminolevulinate synthase added to an organism
   that makes ALA via the C5 pathway (hemA/hemL); forcing the C4 route bypasses
   the essential hem genes and breaks gene-essentiality predictions.
2. FAKE GENE: the attached gene is a WEAK homology hit (pident < 40) whose real
   function is likely a DIFFERENT EC (paralog / annotation artifact).
3. NON-PHYSIOLOGICAL DIRECTION: the reaction is forced to run in a reverse or
   unphysiological direction to rescue growth (e.g. an aminotransferase run
   backwards to synthesize a biosynthetic intermediate, bypassing the committed
   step of the real pathway).
4. WRONG TRANSPORTER: substrate incompatible with the organism's physiology /
   habitat / Gram character.

## INCLUDE when

- The enzyme/transporter is genuinely present (strong BLAST, pident >= 60).
- It is a legitimate step the organism is known to have, filling a true
  annotation gap (enzyme present but unannotated; pathway context supports it).
- Spontaneous reaction (no enzyme needed).
- The reaction carries no gene but is a well-established step of a confirmed
  pathway of this organism.

## GENELESS reactions (reactions marked NO GENE EVIDENCE)

Many real metabolic steps have NO reliable BLAST hit (transmembrane/short
enzymes, species-specific paralogs, or BLAST simply failed to align). These
"geneless" reactions are often REAL, not bypasses.

For a GENELESS reaction, judge by PATHWAY MEMBERSHIP, not by the absence of gene
evidence:
- If its pathway (read from its name/equation) appears in the organism's
  "Predicted Pathways" list in the prompt context, INCLUDE it. This is an
  annotation gap, not a foreign pathway. Examples of near-universal pathways
  that are frequently geneless in models: LPS core synthesis, cardiolipin
  synthesis, thiamine/folate/spermidine biosynthesis, fatty-acid synthesis.
- EXCLUDE a geneless reaction ONLY when its pathway is one the organism is
  KNOWN NOT to use (e.g. a sulfide-dependent thiamine route in an organism that
  uses ThiFSGH; a carboxyspermidine route in an organism that uses SpeD/SpeE).

## Confidence rules

- Use confidence "high" ONLY for a clear foreign pathway or an obvious bypass,
  supported by positive evidence (e.g. the organism is known to use a different
  route). The mere ABSENCE of a BLAST hit is NOT sufficient for "high".
- Anything you are not sure about is "medium" or "low".

Output ONLY valid JSON.
"""


def _load_taxonomy_tree(output_dir: str) -> str:
    """GTDB lineage from the run's taxonomy/taxonomy.tsv (written by code/taxonomy.py)."""
    return taxonomy.lineage_of(output_dir)


def _load_mimedb_profile(output_dir: str) -> str:
    path = os.path.join(output_dir, "prompts", "pathway_verification_system.txt")
    if not os.path.exists(path):
        return ""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    if "[TARGET ORGANISM CONTEXT]" not in text:
        return ""
    start = text.find("[TARGET ORGANISM CONTEXT]")
    end = text.find("\n---", start + 30)
    if end < 0:
        end = text.find("\n[ROLE]", start + 30)
    if end > start:
        return text[start:end].strip()
    return ""


def _load_pathway_inventory(output_dir: str, genome_id: str,
                            max_items: int = 120) -> str:
    """Load the organism's predicted pathway inventory from {genome}-all-Pathways.tbl.

    Returns a compact "ID — Name" list of pathways with Prediction=TRUE, sorted
    by completeness (highest first), bounded to max_items lines.
    """
    candidates = glob.glob(os.path.join(output_dir, "*-all-Pathways.tbl"))
    if not candidates:
        return ""
    path = candidates[0]

    rows = []
    try:
        with open(path, encoding="utf-8", errors="replace") as fh:
            for ln in fh:
                if not ln.strip() or ln.lstrip().startswith("#"):
                    continue
                parts = ln.rstrip("\n").split("\t")
                if len(parts) < 4:
                    continue
                pid, name, prediction, completeness = parts[0], parts[1], parts[2], parts[3]
                # keep only pathways the organism is predicted to have
                if prediction.strip().upper() != "TRUE":
                    continue
                try:
                    comp = float(completeness)
                except (TypeError, ValueError):
                    comp = 0.0
                rows.append((comp, pid, name))
    except OSError:
        return ""

    if not rows:
        return ""

    rows.sort(key=lambda r: -r[0])
    lines = [f"{pid} — {html.unescape(name)}" for _, pid, name in rows[:max_items]]
    return "\n".join(lines)


def _format_reaction(r: dict, seed_db: dict, idx: int, total: int) -> str:
    seed = r.get("seed", "")
    lines = [f"  [{idx}/{total}] Seed: {seed}"]
    db = seed_db.get(seed, {})
    name = r.get("name") or (db.get("name") if db else "") or ""
    if name:
        lines.append(f"      Name:      {name}")
    eq = r.get("equation") or (db.get("definition") if db else "") or ""
    if eq:
        lines.append(f"      Equation:  {eq}")
    ec = r.get("ec")
    if ec:
        lines.append(f"      EC:        {ec}")
    pwy = r.get("pathway")
    pwy_status = r.get("pathway_status")
    if pwy:
        lines.append(f"      Pathway:   {pwy}  (status: {pwy_status or 'none'})")
    bs = r.get("bitscore")
    pident = r.get("pident")
    if isinstance(bs, (int, float)):
        tag = "STRONG" if (bs or 0) >= 200 else ("MODERATE" if (bs or 0) >= 100 else "WEAK")
        lines.append(f"      Bitscore:  {bs:.1f} ({tag}),  pident: {pident if pident is not None else 'NA'}%")
    genes = r.get("genes") or ""
    if genes:
        gene_list = [g for g in str(genes).split("|") if g][:6]
        lines.append(f"      Genes:     {', '.join(gene_list)}")
    if r.get("spontaneous"):
        lines.append("      SPONTANEOUS (no enzyme needed)")
    if r.get("geneless"):
        lines.append("      NO GENE EVIDENCE (judge by PATHWAY MEMBERSHIP against the Predicted Pathways list)")
    return "\n".join(lines)


def _build_prompt(added: list, org_name: str, taxonomy_tree: str,
                  mimedb_profile: str, seed_db: dict,
                  pathway_inventory: str) -> str:
    parts = []
    parts.append("=" * 72)
    parts.append("ORGANISM CONTEXT")
    parts.append("=" * 72)
    parts.append(f"Organism: {org_name}")
    if taxonomy_tree:
        parts.append(f"GTDB-Tk Taxonomy: {taxonomy_tree}")
    if mimedb_profile:
        parts.append(f"\n{mimedb_profile[:2000]}")
    if pathway_inventory:
        parts.append("")
        parts.append(f"Predicted Pathways (from genome annotation, highest-completeness first):")
        parts.append(pathway_inventory[:6000])
    parts.append("")
    parts.append("=" * 72)
    parts.append(f"REACTIONS ADDED BY GAPFILLING ({len(added)})")
    parts.append("=" * 72)
    parts.append("These reactions were added only because FBA required them for biomass.")
    parts.append("")
    for i, r in enumerate(added, 1):
        parts.append(_format_reaction(r, seed_db, i, len(added)))
        parts.append("")
    parts.append("=" * 72)
    parts.append("INSTRUCTIONS")
    parts.append("=" * 72)
    parts.append(
        f"Judge each of the {len(added)} reaction(s) above: INCLUDE if it is a real "
        "metabolic step of this organism, EXCLUDE if it is a wrong-pathway or "
        "bypass leak. For GENELESS reactions, use the Predicted Pathways list to "
        "judge pathway membership. Use 'high' confidence ONLY for clear foreign "
        "pathways/bypasses with positive evidence. Output ONLY valid JSON."
    )
    parts.append("""
===== JSON OUTPUT CONTRACT (strict) =====
Respond with ONE valid JSON object only, no markdown fences, no commentary.
Shape:
{"decisions": [
  {"seed": "<reaction_id>", "decision": "include|exclude",
   "confidence": "high|medium|low", "rationale": "<one-sentence biological reason>"}
]}
Every seed listed above must appear exactly once.
""")
    return "\n".join(parts)


def _parse_decisions(text: str) -> list:
    if not isinstance(text, str) or not text.strip():
        raise ValueError("empty LLM response")
    cleaned = text.strip()
    if cleaned.startswith("```"):
        cleaned = re.sub(r"^```(?:json)?\s*", "", cleaned)
        cleaned = re.sub(r"\s*```$", "", cleaned)
        cleaned = cleaned.strip()
    payload = json.loads(cleaned)
    if isinstance(payload, dict) and "decisions" in payload:
        items = payload["decisions"]
    elif isinstance(payload, list):
        items = payload
    else:
        raise ValueError(f"unexpected JSON shape: {type(payload).__name__}")
    decisions = []
    for item in items:
        if not isinstance(item, dict):
            continue
        seed = item.get("seed")
        decision = item.get("decision", "include")
        if decision not in ("include", "exclude"):
            decision = "include"
        if seed is None:
            continue
        decisions.append({
            "seed": str(seed),
            "decision": decision,
            "confidence": item.get("confidence", "low"),
            "rationale": item.get("rationale") or "",
        })
    if not decisions:
        raise ValueError("no decisions parsed")
    return decisions


def _review(added: list, org_name: str, taxonomy_tree: str,
            mimedb_profile: str, seed_db: dict,
            pathway_inventory: str, max_retries: int = 2,
            batch_size: int = 5) -> list:
    """Ask the LLM in small parallel batches so each JSON response stays well
    under the completion-token limit (10+ reactions per prompt truncated)."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    structured_llm = llm.bind(response_format={"type": "json_object"}, thinking="enabled")
    n_batches = (len(added) + batch_size - 1) // batch_size

    def _one_batch(i):
        batch = added[i:i + batch_size]
        prompt = _build_prompt(batch, org_name, taxonomy_tree, mimedb_profile,
                               seed_db, pathway_inventory)
        print(f"  [review] batch {i // batch_size + 1}/{n_batches} ({len(batch)} reactions)")
        for attempt in range(max_retries + 1):
            try:
                messages = [{"role": "system", "content": GAPFILL_SYSTEM_PROMPT},
                            {"role": "user", "content": prompt}]
                response = structured_llm.invoke(messages)
                decisions = _parse_decisions(response.content)
                expected = {r["seed"] for r in batch}
                got = {d["seed"] for d in decisions}
                missing = expected - got
                for s in missing:
                    decisions.append({
                        "seed": s, "decision": "include", "confidence": "low",
                        "rationale": "LLM did not provide a decision; keep (insufficient evidence to exclude).",
                    })
                return decisions
            except Exception as e:
                print(f"  [RETRY {attempt + 1}/{max_retries}] {e}")
                if attempt < max_retries:
                    time.sleep(2 ** attempt)
                else:
                    return [{
                        "seed": r["seed"], "decision": "include", "confidence": "low",
                        "rationale": f"LLM error: {str(e)[:200]}; keep (insufficient evidence to exclude).",
                    } for r in batch]

    all_decisions = []
    with ThreadPoolExecutor(max_workers=n_batches) as pool:
        futures = [pool.submit(_one_batch, i) for i in range(0, len(added), batch_size)]
        for fut in as_completed(futures):
            all_decisions.extend(fut.result())
    return all_decisions


def _read_exclude(exclude_file: str) -> list:
    if not os.path.exists(exclude_file):
        return []
    seeds = []
    with open(exclude_file, encoding="utf-8") as fh:
        for ln in fh:
            ln = ln.strip()
            if not ln or ln == "seed":
                continue
            seeds.append(ln.split("\t")[0])
    return seeds


def _append_exclude(exclude_file: str, seeds: list) -> int:
    existing = set(_read_exclude(exclude_file))
    new = [s for s in seeds if s not in existing]
    if not new:
        return 0
    with open(exclude_file, "a", encoding="utf-8") as fh:
        for s in new:
            fh.write(s + "\n")
    return len(new)



def _undo_append_exclude(exclude_file: str, seeds: list) -> int:
    """Remove `seeds` from the exclusion list (restore pre-review state)."""
    current = _read_exclude(exclude_file)
    keep = [s for s in current if s not in set(seeds)]
    removed = len(current) - len(keep)
    if removed == 0:
        return 0
    with open(exclude_file, "w", encoding="utf-8") as fh:
        fh.write("seed\n")
        for s in keep:
            fh.write(s + "\n")
    return removed

def _run_gapfill(output_dir: str, genome_id: str, draft_model: str,
                 medium_file: str, bcore: int, exclude_file: str) -> bool:
    """Run official gapseq fill with the review's exclusions.

    The exclusions are applied by filtering gapseq's candidate RDS tables and
    invoking the installed ``gapseq`` CLI; this file never patches or imports
    gapseq's R implementation.
    """
    from gemagent.adapters.gapseq import fill_with_excludes
    print(f"  Re-running official gapseq fill (exclude list: {exclude_file}) ...")
    try:
        return fill_with_excludes(
            output_dir, genome_id, draft_model, medium_file, exclude_file,
            bcore=bcore,
        )
    except Exception as exc:
        print(f"  [ERROR] official gapseq fill failed: {exc}", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output_dir", required=True)
    parser.add_argument("--genome_id", required=True)
    parser.add_argument("--org_name", default="")
    parser.add_argument("--taxonomy_tree", default="")
    parser.add_argument("--draft_model", required=True)
    parser.add_argument("--medium", required=True)
    parser.add_argument("--bcore", type=int, default=100)
    parser.add_argument("--max_iter", type=int, default=3)
    parser.add_argument("--retry-one-by-one", action="store_true",
                        help="If re-gapfilling fails after adding exclusions, retry each "
                             "exclusion individually and keep only the feasible ones "
                             "(instead of reverting the whole batch).")
    args = parser.parse_args()

    exclude_file = os.path.join(args.output_dir, f"{args.genome_id}-gapfill-exclude.tsv")
    added_file = os.path.join(args.output_dir, f"{args.genome_id}-gapfill-added.json")
    decisions_file = os.path.join(args.output_dir, f"{args.genome_id}-gapfill-decisions.json")

    taxonomy_tree = args.taxonomy_tree or _load_taxonomy_tree(args.output_dir)
    mimedb_profile = _load_mimedb_profile(args.output_dir)
    pathway_inventory = _load_pathway_inventory(args.output_dir, args.genome_id)
    seed_db = _load_seed_db()

    all_decisions = {}
    for iteration in range(1, args.max_iter + 1):
        if not os.path.exists(added_file):
            print(f"  [iter {iteration}] no {os.path.basename(added_file)}; stopping.")
            break
        with open(added_file, encoding="utf-8") as fh:
            added = json.load(fh)
        if not isinstance(added, list) or not added:
            print(f"  [iter {iteration}] no gapfill-added reactions; converged.")
            break

        print(f"\n  [iter {iteration}] reviewing {len(added)} gapfill-added reaction(s) ...")
        decisions = _review(added, args.org_name, taxonomy_tree,
                            mimedb_profile, seed_db, pathway_inventory)
        for d in decisions:
            all_decisions[d["seed"]] = d
        with open(decisions_file, "w", encoding="utf-8") as fh:
            json.dump(decisions, fh, indent=2, ensure_ascii=False)

        # OPT-2: only HIGH-confidence EXCLUDE goes to the permanent exclusion list.
        high_excludes = [d["seed"] for d in decisions
                         if d["decision"] == "exclude" and d.get("confidence") == "high"]
        n_exc = sum(1 for d in decisions if d["decision"] == "exclude")
        n_inc = len(decisions) - n_exc
        print(f"  [iter {iteration}] include={n_inc}, exclude={n_exc} "
              f"(high-confidence excludes -> permanent: {len(high_excludes)})")
        for d in decisions:
            if d["decision"] == "exclude":
                tag = "PERM" if d.get("confidence") == "high" else "soft"
                print(f"      EXCLUDE[{tag}] {d['seed']} ({d.get('confidence')}): {d['rationale'][:110]}")

        if not high_excludes:
            print("  No HIGH-confidence bypasses; converged (medium/low excludes are kept).")
            break

        n_added = _append_exclude(exclude_file, high_excludes)
        print(f"  Appended {n_added} HIGH-confidence reaction(s) to exclusion list.")

        if not _run_gapfill(args.output_dir, args.genome_id, args.draft_model,
                            args.medium, args.bcore, exclude_file):
            if args.retry_one_by_one:
                _undo_append_exclude(exclude_file, high_excludes)
                kept = []
                for _seed in high_excludes:
                    _append_exclude(exclude_file, [_seed])
                    if _run_gapfill(args.output_dir, args.genome_id, args.draft_model,
                                    args.medium, args.bcore, exclude_file):
                        kept.append(_seed)
                    else:
                        _undo_append_exclude(exclude_file, [_seed])
                        print(f"      [retry-one-by-one] skip {_seed} (re-gapfill infeasible)",
                              file=sys.stderr)
                print(f"  [retry-one-by-one] kept {len(kept)}/{len(high_excludes)} exclusion(s); "
                      f"final model re-gapfilled.", file=sys.stderr)
                if kept:
                    continue
                break
            n_undone = _undo_append_exclude(exclude_file, high_excludes)
            print(f"  [WARN] re-gapfilling failed after excluding {n_added} reaction(s); "
                  f"reverted {n_undone} exclusion(s). Keeping the pre-review "
                  f"gapfilled model as the accepted output.", file=sys.stderr)
            print("  Decisions are still recorded for reference/manual review.", file=sys.stderr)
            break

    # final consolidated decisions (latest decision per seed wins)
    final_path = os.path.join(args.output_dir, f"{args.genome_id}-gapfill-decisions-final.json")
    with open(final_path, "w", encoding="utf-8") as fh:
        json.dump(list(all_decisions.values()), fh, indent=2, ensure_ascii=False)
    n_exc = sum(1 for d in all_decisions.values() if d["decision"] == "exclude")
    n_exc_high = sum(1 for d in all_decisions.values()
                     if d["decision"] == "exclude" and d.get("confidence") == "high")
    print(f"\n  Gapfill review done: {len(all_decisions)} reactions reviewed, "
          f"{n_exc} excluded ({n_exc_high} high-confidence permanent).")
    return


if __name__ == "__main__":
    main()
