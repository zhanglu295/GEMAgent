#!/usr/bin/env python3
"""Minimal SEED reaction dictionary loader.

Extracted from the retired candidate reviewer (review_candidates.py). It survives
because the post-gapfill bypass review (review_gapfill_added.py) still needs the
SEED reaction metadata (name / definition / EC) to build its prompts.

The table comes from :mod:`paths` rather than a hardcoded ``database/`` path, and
the suppression column is detected from the real header. This matters because the
two SEED tables in the tree differ in exactly that column name
(``code/dat`` -> ``gapseq.status``, ``database/`` -> ``gemagent.status``); pinning
either name broke the loader as soon as it was pointed at the other tree.
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import paths  # noqa: E402


def _seed_reaction_paths() -> tuple[str, str]:
    """``(tsv_path, status_column)`` for the SEED reaction table."""
    return paths.seed_reactions_path()


#: Kept for backwards compatibility with anything importing the old constant.
SEED_RXN_TSV, _STATUS_COL = _seed_reaction_paths()


def _load_seed_db() -> dict:
    """SEED reaction ID -> {name, definition, ec}."""
    path, status_col = _seed_reaction_paths()
    if not os.path.exists(path):
        print(f"  [WARN] SEED reaction DB not found: {path}")
        return {}

    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if status_col:
        # Keep rows the pipeline has not suppressed. `removed` is the only
        # suppression value the curated tree uses; anything else stays.
        df = df[df[status_col].str.strip().str.lower() != "removed"]
    else:
        print(f"  [WARN] no status column in {os.path.basename(path)}; "
              f"keeping all {len(df)} rows")

    # NOTE: ec_numbers is empty in BOTH trees (verified 2026-09-19: 0 of 34936
    # rows carry a value), so `ec` is almost always "". Callers should treat it
    # as best-effort; the KO definitions embedded in the GPR-review prompts are
    # the reliable EC source.
    cols = [c for c in ("name", "definition", "ec_numbers") if c in df.columns]
    recs = df.set_index("id")[cols].to_dict("index")
    db = {
        rid: {
            "name": rec.get("name", ""),
            "definition": rec.get("definition", ""),
            "ec": rec.get("ec_numbers", ""),
        }
        for rid, rec in recs.items()
    }
    print(f"  Loaded {len(db)} SEED reactions (from {os.path.basename(os.path.dirname(path))})")
    return db


load_seed_db = _load_seed_db