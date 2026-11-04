#!/usr/bin/env python3
"""us_reader.py — US review v3, stage 1: LLM phenotype reader.

Replaces the rule-based candidate scanner (us_spans_v23.py).  Mechanism is the
deepresearch one (sequential windows over the full text + verbatim-quote
verification), with a phenotype schema:

    consumed | not_consumed | secreted | not_secreted | required
      -> (polarity, action):  consumed/required -> (+, uptake)
                              not_consumed      -> (-, uptake)
                              secreted          -> (+, secretion)
                              not_secreted      -> (-, secretion)

Output: one row per statement, TSV, ready for the review controller's gate.

usage:
  python us_reader.py --papers <papers.json> --org "Bifidobacterium adolescentis" \
      [--targets targets.txt] --out claims_llm.tsv [--win 24000] [--limit 3]
"""
from __future__ import annotations

import argparse, json, os, re, sys, tempfile, time
from pathlib import Path
from typing import Dict, List

import pandas as pd

AGENT = os.environ.get("GEMAGENT_AGENT_DIR", "").strip()
DEEP_RESEARCH_DIR = os.environ.get("GEMAGENT_DEEP_RESEARCH_DIR", "").strip()
if not DEEP_RESEARCH_DIR and AGENT:
    DEEP_RESEARCH_DIR = str(Path(AGENT) / "agent_deepresearch")
for _p in (DEEP_RESEARCH_DIR, AGENT):
    if _p and _p not in sys.path:
        sys.path.insert(0, _p)

try:
    import extract as X                                  # noqa: E402
    from dr_evidence_utils import (                      # noqa: E402
        clean_org_name, organism_tokens, split_sentences, verify_verbatim)
except ImportError as exc:
    raise ImportError(
        "Deep Research reader helpers were not found. Set "
        "GEMAGENT_DEEP_RESEARCH_DIR to the directory containing extract.py and "
        "dr_evidence_utils.py."
    ) from exc

BALANCE_MARKER = Path(os.environ.get(
    "GEMAGENT_BALANCE_MARKER",
    str(Path(tempfile.gettempdir()) / "gemagent_reader" / ".BALANCE_EXHAUSTED"),
))

CAT2PA = {
    "consumed":     ("positive", "uptake"),
    "required":     ("positive", "uptake"),
    "not_consumed": ("negative", "uptake"),
    "secreted":     ("positive", "secretion"),
    "not_secreted": ("negative", "secretion"),
}

SYSTEM_TMPL = """You are a metabolic-model curator reading a paper to find PHENOTYPE statements about ONE organism.

[TARGET ORGANISM]
{org}
(also printed as: {aliases})

[TASK]
From the supplied text window, extract every statement that ties a metabolite to the TARGET organism, using exactly one of these five category words:
  "consumed"      the organism takes it up / grows on it / uses it as a substrate
  "not_consumed"  the organism explicitly cannot use it (no growth, negative for utilization, unable to use X)
  "secreted"      the organism produces and releases it (fermentation product, end product, excretion)
  "not_secreted"  the organism explicitly does not produce/release it (negative for production)
  "required"      it must be supplied (growth requirement, auxotrophy, supplement) -- i.e. the organism has to take it up

[RULES]
- The quote MUST be copied verbatim from the supplied text (exact string, one sentence or clause). Never paraphrase; never join two sentences.
- Only statements about the TARGET organism count. Results obtained in a different organism are NOT evidence for the target.
- subject_level: "strain" (this strain / an isolate identified as it), "species" (the species, or a named strain of it), "genus" (the genus / another species of it), "other" (unrelated organism, medium recipe, engineered strain).
- evidence_form: "phenotype_test" (a direct test or observation of the organism: grows / does not grow / utilized / negative for ...), "inference" (deduced, e.g. from genome or modelling), "other".
- A compound that is merely mentioned (buffers, medium components, analytical reagents) must NOT be reported unless the text asserts a phenotype of the target organism.
- Report one entry per compound statement.

[OUTPUT]
Return ONLY JSON:
{{"phenotypes": [{{"compound": "Acetate",
  "category": "consumed|not_consumed|secreted|not_secreted|required",
  "quote": "<verbatim>",
  "organism_mentioned": "<organism string as printed, or empty>",
  "subject_level": "strain|species|genus|other",
  "evidence_form": "phenotype_test|inference|other",
  "confidence": "high|medium|low"}}]}}
Return {{"phenotypes": []}} when the window holds no target-specific evidence."""

TARGETS_BLOCK = """

[COMPOUNDS SPECIFICALLY UNDER INVESTIGATION]
{targets}

These are the phenotypes this run exists to resolve. If the text contains a verbatim statement that ties ANY compound from that list to the TARGET organism -- consumed, not consumed, secreted, not secreted, or required -- report it, even if it is only one clause of a longer sentence. Report every other metabolite statement you find as well: the list adds to the task, it does not restrict it."""


def load_targets(path: str) -> List[str]:
    if not path or not os.path.exists(path):
        return []
    out = []
    for line in open(path, encoding="utf-8"):
        t = line.strip()
        if t and not t.startswith("#"):
            out.append(t)
    return out


def build_system(org: str, targets: List[str]) -> str:
    toks = organism_tokens(org)
    aliases = ", ".join(sorted(toks.get("genus_aliases") or ())) or "(none)"
    base = SYSTEM_TMPL.replace("{org}", org).replace("{aliases}", aliases)
    if targets:
        base += TARGETS_BLOCK.replace("{targets}", ", ".join(targets))
    return base


def build_user(org: str, paper: dict, window: str, part: int, parts: int,
               targets: List[str]) -> str:
    toks = organism_tokens(org)
    strains = ", ".join(toks.get("strain_tokens") or ()) or "(none)"
    accs = ", ".join(toks.get("accessions") or ()) or "(none)"
    wanted = ("\nCompounds this run is specifically looking for: " + ", ".join(targets)
              + " -- report any of them the text ties to the organism.\n") if targets else ""
    return f"""Target organism: {org}
Accepted names: {clean_org_name(org)}
Strain tokens: {strains}
Assembly accessions: {accs}
{wanted}
Paper: {paper.get('title')} ({paper.get('year') or 'n.d.'})
Text window {part}/{parts}:

{window}

Extract every phenotype statement about the TARGET organism (consumed / not_consumed / secreted / not_secreted / required). Quote verbatim. Return strict JSON only."""


def attr_level(org: str, quote: str, organism_mentioned: str) -> str:
    toks = organism_tokens(org)
    ql = (quote or "").lower()
    om = (organism_mentioned or "").lower()
    for t in toks.get("strain_tokens") or []:
        if t and (t.lower() in ql or t.lower() in om):
            return "strain"
    sp = clean_org_name(org).lower()
    if sp and (sp in ql or sp in om):
        return "species"
    for a in toks.get("genus_aliases") or []:
        if a and (a.lower() in ql or a.lower() in om):
            return "species"
    genus = (org.split()[0] if org else "").lower()
    if genus and genus in ql:
        return "genus"
    return "none"


def read_paper(paper: dict, org: str, llm, system_prompt: str, targets: List[str],
               win_chars: int) -> List[dict]:
    text = paper.get("text") or ""
    if not text.strip():
        return []
    parts = X.windows_of(text, win_chars)
    out: List[dict] = []
    for idx, window in enumerate(parts, 1):
        prompt = build_user(org, paper, window, idx, len(parts), targets)
        raw = ""
        for attempt in range(3):
            try:
                resp = llm.bind(response_format={"type": "json_object"}).invoke(
                    [{"role": "system", "content": system_prompt},
                     {"role": "user", "content": prompt}], max_tokens=8192)
                raw = getattr(resp, "content", "") or ""
            except Exception as exc:                                  # noqa: BLE001
                print(f"      LLM error (try {attempt+1}): {exc}")
                if "Insufficient Balance" in str(exc) or "Error code: 402" in str(exc):
                    BALANCE_MARKER.parent.mkdir(parents=True, exist_ok=True)
                    BALANCE_MARKER.write_text(str(exc), encoding="utf-8")
                    raise RuntimeError("LLM balance exhausted (402): aborting run") from exc
                raw = ""
            if raw.strip():
                break
            time.sleep(1 + attempt)
        if not raw.strip():
            continue
        payload = X._parse_json(raw)
        items = payload.get("phenotypes") or payload.get("metabolites") or []
        kept = dropped = 0
        for it in items:
            if not isinstance(it, dict):
                continue
            name = (it.get("compound") or it.get("name") or "").strip()
            cat = str(it.get("category") or "").strip().lower()
            quote = (it.get("quote") or "").strip()
            if not name or not quote or cat not in CAT2PA:
                dropped += 1
                continue
            ok = verify_verbatim(quote, window) or verify_verbatim(quote, text)
            if not ok:
                dropped += 1
                continue
            pol, act = CAT2PA[cat]
            out.append(dict(
                compound=name, category=cat, polarity=pol, action=act,
                required=("required" == cat),
                subject_level=(it.get("subject_level") or "").strip().lower(),
                evidence_form=(it.get("evidence_form") or "").strip().lower(),
                confidence=(it.get("confidence") or "").strip().lower(),
                attr=attr_level(org, quote, it.get("organism_mentioned") or ""),
                verbatim=1,
                quote=quote,
                doi=(paper.get("doi") or ""),
                paper=(paper.get("title") or "")[:120],
                window=idx,
            ))
            kept += 1
        print(f"    win {idx}/{len(parts)}: kept {kept}, dropped {dropped}")
    # dedupe
    seen, uniq = set(), []
    for r in out:
        key = (r["compound"].lower(), r["category"], r["quote"])
        if key in seen:
            continue
        seen.add(key)
        uniq.append(r)
    return uniq


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--papers", required=True)
    ap.add_argument("--org", required=True)
    ap.add_argument("--targets", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--win", type=int, default=24000)
    ap.add_argument("--limit", type=int, default=0, help="only first N papers (smoke test)")
    a = ap.parse_args()

    papers = json.load(open(a.papers, encoding="utf-8"))
    if a.limit:
        papers = papers[:a.limit]
    targets = load_targets(a.targets)
    system_prompt = build_system(a.org, targets)
    llm = X._load_llm()
    if llm is None:
        print("[us_reader] no LLM available"); return 1

    print(f"[us_reader] org={a.org} | papers={len(papers)} | targets={len(targets)} | win={a.win}")
    rows: List[dict] = []
    t0 = time.time()
    for i, paper in enumerate(papers, 1):
        ttl = (paper.get("title") or "")[:70]
        print(f"  [{i}/{len(papers)}] {ttl}")
        rows += read_paper(paper, a.org, llm, system_prompt, targets, a.win)
    df = pd.DataFrame(rows)
    df.insert(0, "strain", a.org)
    df.to_csv(a.out, sep="\t", index=False)
    print(f"[us_reader] wrote {a.out}: {len(df)} claims | {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
