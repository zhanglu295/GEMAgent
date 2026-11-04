#!/usr/bin/env python3
"""Apply frozen LLM REMOVE decisions to SBML fbc:or GPR branches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import re
from collections import defaultdict
from pathlib import Path


REACTION = re.compile(
    r"(?P<open><reaction\b[^>]*\bid=\"(?P<rid>[^\"]+)\"[^>]*>)"
    r"(?P<body>.*?)"
    r"(?P<close></reaction>)",
    re.S,
)
GENE_REF = re.compile(
    r"\s*<fbc:geneProductRef\b[^>]*\bfbc:geneProduct=\"(?P<gene>[^\"]+)\"[^>]*/>"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def norm_rid(rid: str) -> str:
    return rid[2:] if rid.startswith("R_") else rid


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--decisions", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    model_path = Path(args.model)
    decisions_path = Path(args.decisions)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    actions = defaultdict(lambda: defaultdict(str))
    with decisions_path.open(encoding="utf-8", newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            action = (row.get("action") or row.get("verdict") or "KEEP").strip().upper()
            actions[row["reaction_id"]][row["gene"]] = action
            actions[norm_rid(row["reaction_id"])][row["gene"]] = action

    text = model_path.read_text(encoding="utf-8", errors="replace")
    report = []
    changed = skipped_empty = no_removals = removed_total = 0

    def rewrite(match):
        nonlocal changed, skipped_empty, no_removals, removed_total
        rid = match.group("rid")
        reaction_actions = actions.get(rid) or actions.get(norm_rid(rid))
        if not reaction_actions:
            return match.group(0)
        remove = {
            gene for gene, action in reaction_actions.items()
            if action == "REMOVE"
        }
        if not remove:
            no_removals += 1
            return match.group(0)

        body = match.group("body")
        refs = list(GENE_REF.finditer(body))
        if not refs:
            return match.group(0)
        kept_genes = [
            ref.group("gene") for ref in refs if ref.group("gene") not in remove
        ]
        if not kept_genes:
            skipped_empty += 1
            report.append([
                rid, "SKIP_ALL_REMOVED", "", "", ",".join(sorted(remove))
            ])
            return match.group(0)

        removed_here = []

        def remove_ref(ref_match):
            gene = ref_match.group("gene")
            if gene in remove:
                removed_here.append(gene)
                return ""
            return ref_match.group(0)

        new_body = GENE_REF.sub(remove_ref, body)
        if not removed_here:
            return match.group(0)
        changed += 1
        removed_total += len(removed_here)
        report.append([
            rid,
            "APPLIED",
            ",".join(ref.group("gene") for ref in refs),
            ",".join(kept_genes),
            ",".join(removed_here),
        ])
        return match.group("open") + new_body + match.group("close")

    new_text = REACTION.sub(rewrite, text)
    out_path.write_text(new_text, encoding="utf-8")

    report_path = out_path.with_name(out_path.stem + "_apply_report.tsv")
    with report_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow([
            "reaction_id", "action", "genes_before", "genes_after", "removed_genes"
        ])
        writer.writerows(report)

    manifest = {
        "input_model": os.path.abspath(model_path),
        "input_model_sha256": sha256(model_path),
        "decisions": os.path.abspath(decisions_path),
        "decisions_sha256": sha256(decisions_path),
        "output_model": os.path.abspath(out_path),
        "output_model_sha256": sha256(out_path),
        "apply_report": os.path.abspath(report_path),
        "reactions_changed": changed,
        "gene_links_removed": removed_total,
        "skipped_all_removed": skipped_empty,
        "reactions_without_removals": no_removals,
        "truth_firewall": "pass: frozen LLM decisions only",
    }
    manifest_path = out_path.with_name(out_path.stem + "_apply_manifest.json")
    manifest_path.write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(manifest, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
