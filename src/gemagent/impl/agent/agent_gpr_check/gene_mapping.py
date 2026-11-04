#!/usr/bin/env python3
"""Deterministic gene-ID resolution for GPR evidence.

A model gene may be encoded as:
  - a KOfam/GFF gene ID:       NC_000913.3_1710
  - that ID plus a local ORF:  NC_000913_3_1710_1_753
  - genomic coordinates:        NC_000913_3_4449652_4449125
  - an external gene with a separate coordinate table (CarveMe/ModelSEED).

Resolution is deliberately fail-closed: aliases must identify one KOfam gene.
No generic "last number" fallback is allowed because the final number can be a
local fragment end rather than the gene ordinal.
"""
from __future__ import annotations

import bisect
import gzip
import re
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

LOCAL_COORD_RE = re.compile(r"^(.+)_(-?\d+)[:_](-?\d+)$")
TRAILING_NUM_RE = re.compile(r"(\d+)$")


def norm_alias(value: object) -> str:
    s = str(value).strip()
    if s.startswith("G_"):
        s = s[2:]
    return re.sub(r"[^A-Za-z0-9_]", "_", s)


def aliases(value: object) -> Set[str]:
    s = str(value).strip()
    return {norm_alias(s), norm_alias(s[2:] if s.startswith("G_") else s)}


def last_num(value: object) -> Optional[str]:
    m = TRAILING_NUM_RE.search(str(value).strip())
    return m.group(1) if m else None


def open_text(path: str | Path):
    p = Path(path)
    if p.suffix == ".gz":
        return gzip.open(p, "rt", encoding="utf-8", errors="replace")
    return p.open(encoding="utf-8", errors="replace")


def load_kofam_index(path: str | Path) -> dict:
    by_id: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    by_alias: Dict[str, Set[str]] = defaultdict(set)
    by_num: Dict[str, Set[str]] = defaultdict(set)
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split(None, 6)
            if len(parts) < 7 or parts[0] != "*":
                continue
            gene_id = parts[1].strip()
            if not gene_id:
                continue
            by_id[gene_id].append((parts[2].strip(), parts[6].strip().strip('"')))
            for alias in aliases(gene_id):
                by_alias[alias].add(gene_id)
            gn = last_num(gene_id)
            if gn:
                by_num[gn].add(gene_id)
    return {
        "by_id": dict(by_id),
        "by_alias": {k: frozenset(v) for k, v in by_alias.items()},
        "by_num": {k: frozenset(v) for k, v in by_num.items()},
    }


def load_gff_records(path: str | Path) -> List[dict]:
    records: List[dict] = []
    if not path:
        return records
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9 or c[2] != "CDS":
                continue
            m = re.search(r"(?:^|;)ID=([^;]+)", c[8])
            if not m:
                continue
            try:
                start, end = int(c[3]), int(c[4])
            except ValueError:
                continue
            records.append({
                "gene_id": m.group(1),
                "contig": c[0],
                "start": min(start, end),
                "end": max(start, end),
            })
    return records


def load_faa_records(path: str | Path) -> List[dict]:
    records: List[dict] = []
    if not path:
        return records
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            parts = line[1:].rstrip("\n").split("#")
            if len(parts) < 3:
                continue
            gene_id = parts[0].strip()
            try:
                start, end = int(parts[1]), int(parts[2])
            except ValueError:
                continue
            records.append({
                "gene_id": gene_id,
                "contig": gene_id.rsplit("_", 1)[0],
                "start": min(start, end),
                "end": max(start, end),
            })
    return records


def _unique_alias_match(identifier: object, kofam: Mapping) -> Optional[str]:
    hits = set()
    for a in aliases(identifier):
        hits.update(kofam["by_alias"].get(a, ()))
    return next(iter(hits)) if len(hits) == 1 else None


def _record_alias_match(identifier: object, records: Sequence[dict]) -> Optional[dict]:
    alias_set = aliases(identifier)
    hits = [r for r in records if aliases(r["gene_id"]) & alias_set]
    if len(hits) == 1:
        return hits[0]
    return None


def _resolve_record_to_kofam(record: Mapping, kofam: Mapping) -> Optional[str]:
    direct = _unique_alias_match(record["gene_id"], kofam)
    if direct:
        return direct
    gn = last_num(record["gene_id"])
    if gn:
        hits = kofam["by_num"].get(gn, ())
        if len(hits) == 1:
            return next(iter(hits))
    return None


def _interval_match(contig: str, start: int, end: int, records: Sequence[dict]) -> Optional[dict]:
    if not records:
        return None
    nc = norm_alias(contig)
    candidates = []
    for r in records:
        ov = max(0, min(end, r["end"]) - max(start, r["start"]))
        if ov <= 0:
            continue
        same_contig = norm_alias(r["contig"]) == nc
        candidates.append((int(same_contig), ov, r))
    if not candidates:
        return None
    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_same, best_ov, best = candidates[0]
    tied = [r for same, ov, r in candidates if same == best_same and ov == best_ov]
    if len({r["gene_id"] for r in tied}) != 1:
        return None
    return best


def parse_literal_coord(identifier: object) -> Optional[Tuple[str, int, int]]:
    s = str(identifier).strip()
    if s.startswith("G_"):
        s = s[2:]
    m = LOCAL_COORD_RE.match(s)
    if not m:
        return None
    contig, a, b = m.group(1), int(m.group(2)), int(m.group(3))
    return contig, min(a, b), max(a, b)


def load_gene_pos_table(path: str | Path | None) -> Dict[str, Tuple[str, int, int]]:
    out: Dict[str, Tuple[str, int, int]] = {}
    if not path:
        return out
    with open_text(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 4:
                c = line.split()
            if len(c) < 4:
                continue
            g, contig, s, e = c[0], c[1], c[2], c[3]
            if str(g).strip().lower() in {"gene", "gene_id", "id"}:
                continue
            try:
                s, e = int(float(s)), int(float(e))
            except ValueError:
                continue
            val = (str(contig).strip(), min(s, e), max(s, e))
            out.setdefault(str(g).strip(), val)
            for a in aliases(g):
                out.setdefault(a, val)
    return out


def resolve_gene_ids(
    model_gene_ids: Iterable[str],
    kofam: Mapping,
    annotation_records: Sequence[dict],
    gene_pos: Optional[Mapping[str, Tuple[str, int, int]]] = None,
) -> Tuple[Dict[str, str], Dict[str, str], Dict[str, str]]:
    """Return model_gene -> KOfam gene, method, and unresolved reason."""
    pos = gene_pos or {}
    mapping: Dict[str, str] = {}
    methods: Dict[str, str] = {}
    unresolved: Dict[str, str] = {}

    for raw in model_gene_ids:
        gid = str(raw)
        direct = _unique_alias_match(gid, kofam)
        if direct:
            mapping[gid], methods[gid] = direct, "kofam_alias"
            continue

        local = LOCAL_COORD_RE.match(gid[2:] if gid.startswith("G_") else gid)
        if local:
            prefix = local.group(1)
            prefix_hit = _unique_alias_match(prefix, kofam)
            if prefix_hit:
                mapping[gid], methods[gid] = prefix_hit, "local_suffix_prefix"
                continue
            record = _record_alias_match(prefix, annotation_records)
            if record is not None:
                hit = _resolve_record_to_kofam(record, kofam)
                if hit:
                    mapping[gid], methods[gid] = hit, "local_suffix_gff"
                    continue

        record = _record_alias_match(gid, annotation_records)
        if record is not None:
            hit = _resolve_record_to_kofam(record, kofam)
            if hit:
                mapping[gid], methods[gid] = hit, "gff_alias"
                continue

        coordinate = None
        for key in (gid, norm_alias(gid)):
            if key in pos:
                coordinate = pos[key]
                source = "gene_pos"
                break
        if coordinate is None:
            parsed = parse_literal_coord(gid)
            if parsed is not None:
                coordinate, source = parsed, "literal_coord"

        if coordinate is not None:
            record = _interval_match(coordinate[0], coordinate[1], coordinate[2], annotation_records)
            if record is not None:
                hit = _resolve_record_to_kofam(record, kofam)
                if hit:
                    mapping[gid], methods[gid] = hit, source
                    continue
            unresolved[gid] = "coordinate_unmatched"
            continue

        unresolved[gid] = "no_unique_evidence_id"

    return mapping, methods, unresolved


def gene_assertions(model_gene_id: str, kofam: Mapping, mapping: Mapping[str, str]) -> List[Tuple[str, str]]:
    key = mapping.get(str(model_gene_id))
    return list(kofam["by_id"].get(key, [])) if key else []
