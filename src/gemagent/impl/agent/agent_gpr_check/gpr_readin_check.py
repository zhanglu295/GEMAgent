#!/usr/bin/env python3
"""Assert that a model actually carries GPRs before and after a GPR review.

The failure this guards against is silent: a CarveMe-style model stores GPRs in
``<notes>GENE_ASSOCIATION: ...`` rather than ``fbc:geneProductAssociation``. A
reader that only understands the fbc form returns zero GPRs, edits nothing, and
still exits 0 -- the review looks successful while doing nothing.

The count is deliberately structural (how many reactions carry an association),
not semantic, so it works for any SBML dialect.

Usage:
    python gpr_readin_check.py --model checked.xml --json post.json --expect-gpr 2170
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

FBC_ASSOC = re.compile(r"<(?:fbc:)?geneProductAssociation\b")
FBC_PRODUCT = re.compile(r'<(?:fbc:)?geneProduct\b[^>]*\bid="([^"]+)"')
NOTES_ASSOC = re.compile(r"GENE_ASSOCIATION:\s*([^<]+)")


def inspect(path: Path) -> dict:
    text = path.read_text(encoding="utf-8", errors="replace")
    fbc = len(FBC_ASSOC.findall(text))
    products = set(FBC_PRODUCT.findall(text))
    notes = len(NOTES_ASSOC.findall(text))
    if fbc:
        mode, gprs = "fbc", fbc
    elif notes:
        mode, gprs = "notes", notes
    else:
        mode, gprs = "none", 0
    return {
        "model": str(path),
        "mode": mode,
        "gpr_reactions": gprs,
        "gene_products": len(products),
        "notes_associations": notes,
        "fbc_associations": fbc,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True)
    ap.add_argument("--json", default="", help="write the summary here")
    ap.add_argument("--expect-gpr", type=int, default=-1,
                    help="fail when the GPR count differs from this (>=0)")
    args = ap.parse_args()

    path = Path(args.model)
    if not path.is_file():
        print(f"gpr_readin_check: missing model {path}", file=sys.stderr)
        return 2
    report = inspect(path)
    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"gpr_readin_check: mode={report['mode']} "
          f"gpr_reactions={report['gpr_reactions']} "
          f"gene_products={report['gene_products']}")

    if report["gpr_reactions"] == 0:
        print("gpr_readin_check: model carries no GPRs; the reader would edit "
              "nothing (silent no-op).", file=sys.stderr)
        return 1
    if args.expect_gpr >= 0 and report["gpr_reactions"] != args.expect_gpr:
        print(f"gpr_readin_check: expected {args.expect_gpr} GPR reactions, "
              f"found {report['gpr_reactions']}.", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
