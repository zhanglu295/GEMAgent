#!/usr/bin/env python3
"""Repair invalid SBML formal charges before reading a model with sybilSBML.

COBRA writers can emit ``charge="-2147483648"`` for metabolites with unknown
charge. That value is outside R's integer range and can become ``NA_integer_``
inside sybilSBML, interrupting model loading. This utility replaces invalid
charges with the SEED value when one is available and otherwise writes zero.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from pathlib import Path

SENTINEL = "-2147483648"
INT32_MIN, INT32_MAX = -2147483647, 2147483647
CHARGE_PATTERN = re.compile(r'(\s)((?:fbc:)?charge)="([^"]*)"')


def default_seed_table() -> Path:
    """Resolve the curated SEED metabolite table from the environment."""
    explicit = os.environ.get("GEMAGENT_SEED_METABOLITES", "").strip()
    if explicit:
        return Path(explicit).expanduser()
    db_dir = os.environ.get("GEMAGENT_DB_DIR", "").strip()
    if db_dir:
        return Path(db_dir).expanduser() / "metabolite" / "seed_metabolites_edited.tsv"
    raise FileNotFoundError(
        "Set GEMAGENT_SEED_METABOLITES or GEMAGENT_DB_DIR before running fix_charges.py"
    )


def load_charge(path: Path) -> dict[str, str]:
    """Load valid integer charges keyed by SEED metabolite identifier."""
    with path.open(encoding="utf-8", errors="replace") as handle:
        rows = [line.rstrip("\n") for line in handle
                if line.strip() and not line.startswith("#")]
    charges: dict[str, str] = {}
    for row in csv.DictReader(rows, delimiter="\t"):
        value = (row.get("charge") or "").strip()
        if not re.fullmatch(r"-?\d+", value):
            continue
        numeric = int(value)
        if INT32_MIN <= numeric <= INT32_MAX:
            charges[row["id"].strip()] = value
    return charges


def fix(path: Path, charge_of: dict[str, str]) -> dict[str, int]:
    """Repair invalid charges in one SBML file and return change counts."""
    xml = path.read_text(encoding="utf-8")
    stats = {"db": 0, "zero": 0, "other": 0}

    def replace(match: re.Match[str]) -> str:
        value = match.group(3)
        valid = re.fullmatch(r"-?\d+", value) and INT32_MIN <= int(value) <= INT32_MAX
        if valid:
            return match.group(0)

        species_start = xml.rfind("<species", 0, match.start())
        species_tag = xml[species_start:match.start()] if species_start >= 0 else ""
        metabolite_match = re.search(r'\bid="([^"]+)"', species_tag)
        metabolite_id = re.sub(r"\[.*\]$", "", metabolite_match.group(1)) if metabolite_match else ""
        if metabolite_id in charge_of:
            stats["db"] += 1
            return f'{match.group(1)}{match.group(2)}="{charge_of[metabolite_id]}"'

        stats["zero" if value == SENTINEL else "other"] += 1
        return f'{match.group(1)}{match.group(2)}="0"'

    repaired = CHARGE_PATTERN.sub(replace, xml)
    if repaired != xml:
        backup = path.with_name(path.name + ".bak_charge")
        if not backup.exists():
            shutil.copy2(path, backup)
        path.write_text(repaired, encoding="utf-8")
    return stats


def main() -> int:
    parser = argparse.ArgumentParser(description="Repair invalid formal charges in SBML files.")
    parser.add_argument("models", nargs="+", help="SBML files to repair in place")
    parser.add_argument(
        "--seed-table",
        default="",
        help="SEED metabolite table; defaults to GEMAGENT_SEED_METABOLITES or GEMAGENT_DB_DIR",
    )
    args = parser.parse_args()

    seed_table = Path(args.seed_table).expanduser() if args.seed_table else default_seed_table()
    if not seed_table.is_file():
        parser.error(f"SEED metabolite table not found: {seed_table}")
    charge_of = load_charge(seed_table)
    print(f"Valid SEED charge records: {len(charge_of)}")

    for raw_path in args.models:
        path = Path(raw_path).expanduser()
        if not path.is_file():
            print(f"  skipped (missing): {path}")
            continue
        stats = fix(path, charge_of)
        print(
            f"  {path}\n"
            f"    repaired_from_seed={stats['db']} "
            f"sentinel_to_zero={stats['zero']} other_to_zero={stats['other']}"
        )
    print("FIX_CHARGE_DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
