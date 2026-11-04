#!/usr/bin/env python3
"""Apply OR2 REMOVE-only decisions to notes-style SBML (CarveMe et al.).

The reaction GPR lives in ``<notes>`` as ``GENE_ASSOCIATION: a or b``. This
applier rewrites only that value and leaves every other byte of the XML
untouched, which avoids cobra round-trip damage on CarveMe models.
"""
import argparse
import csv
import html
import os
import re
from collections import defaultdict
from xml.sax.saxutils import escape

try:
    from .gpr_rewrite import apply_actions
except ImportError:  # script mode
    from gpr_rewrite import apply_actions

REACTION = re.compile(r"<reaction\b([^>]*)>(.*?)</reaction>", re.S)
ID_ATTR = re.compile(r'\bid="([^"]+)"')
GPR = re.compile(r"(GENE_ASSOCIATION:\s*)([^<\r\n]*)", re.I)


def norm_rid(rid):
    """CarveMe SBML prefixes reaction ids with R_; prompts use the bare id."""
    return re.sub(r"^R_", "", str(rid))


def load_decisions(path):
    out = defaultdict(dict)
    with open(path, encoding="utf-8") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            action = (row.get("action") or row.get("verdict") or "KEEP").strip().upper()
            if action != "REMOVE":
                action = "KEEP"
            rid = row["reaction_id"]
            out[rid][row["gene"]] = action
            out.setdefault(norm_rid(rid), {})[row["gene"]] = action
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--model", required=True)
    ap.add_argument("--decisions", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    rxn2actions = load_decisions(args.decisions)
    with open(args.model, encoding="utf-8", errors="replace") as f:
        text = f.read()

    report = []
    n_changed = n_skip_empty = n_no_gpr = n_not_or = 0
    removed_total = added_total = 0

    def rewrite(match):
        nonlocal n_changed, n_skip_empty, n_no_gpr, n_not_or
        nonlocal removed_total, added_total
        attrs, body = match.group(1), match.group(2)
        rid_m = ID_ATTR.search(attrs)
        if not rid_m:
            return match.group(0)
        rid = rid_m.group(1)
        actions = rxn2actions.get(rid) or rxn2actions.get(norm_rid(rid))
        if not actions:
            return match.group(0)

        gpr_m = GPR.search(body)
        if not gpr_m:
            n_no_gpr += 1
            return match.group(0)

        raw = html.unescape(gpr_m.group(2)).strip()
        new, removed, added, skipped = apply_actions(raw, actions)
        if skipped:
            n_skip_empty += 1
            report.append([rid, "SKIP_ALL_REMOVED", raw, raw, "", ""])
            return match.group(0)
        if new == raw:
            n_not_or += 1
            return match.group(0)

        new_body = body[:gpr_m.start(2)] + escape(new) + body[gpr_m.end(2):]
        n_changed += 1
        removed_total += len(removed)
        added_total += len(added)
        report.append([rid, "APPLIED", raw, new, ",".join(removed), ",".join(added)])
        return f"<reaction{attrs}>{new_body}</reaction>"

    new_text = REACTION.sub(rewrite, text)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w", encoding="utf-8") as f:
        f.write(new_text)

    print(f"saved: {args.out}")
    print(f"  reactions changed: {n_changed}")
    print(f"  genes removed: {removed_total}")
    print(f"  genes added: {added_total}")
    print(f"  skipped (all removed): {n_skip_empty}")
    print(f"  no GENE_ASSOCIATION: {n_no_gpr}")
    print(f"  unchanged (non-OR or no-op): {n_not_or}")

    rep_path = os.path.splitext(args.out)[0] + "_apply_report.tsv"
    with open(rep_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["reaction_id", "action", "gpr_before", "gpr_after",
                    "removed_genes", "added_genes"])
        w.writerows(report)
    print(f"  report: {rep_path}")


if __name__ == "__main__":
    main()
