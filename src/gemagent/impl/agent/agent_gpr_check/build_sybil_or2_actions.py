#!/usr/bin/env python3
"""Build sybil-applier actions from frozen OR2 LLM decisions."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--audit_reactions", required=True)
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    reactions = {
        row["reaction_id"]: row
        for row in read_tsv(Path(args.audit_reactions))
    }
    removes = defaultdict(set)
    for row in read_tsv(Path(args.decisions)):
        if row["verdict"].strip().upper() == "REMOVE":
            removes[row["reaction_id"]].add(row["gene"])

    out_rows = []
    for rid, removed in sorted(removes.items()):
        row = reactions.get(rid)
        if row is None:
            continue
        genes = [gene for gene in row["genes"].split(",") if gene]
        kept = [gene for gene in genes if gene not in removed]
        if not kept:
            raise SystemExit(f"refusing to empty GPR: {rid}")
        if len(kept) == len(genes):
            continue
        out_rows.append({
            "reaction_id": rid,
            "gpr_after": " or ".join(kept),
            "genes_before": ",".join(genes),
            "genes_after": ",".join(kept),
            "removed_genes": ",".join(gene for gene in genes if gene in removed),
            "n_before": len(genes),
            "n_after": len(kept),
        })

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(out_rows[0]))
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"wrote {len(out_rows)} actions to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
