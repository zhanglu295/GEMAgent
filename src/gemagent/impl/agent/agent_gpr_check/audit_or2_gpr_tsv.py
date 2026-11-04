#!/usr/bin/env python3
"""Audit pure-OR GPR rows from a sybil model dump with corrected gene IDs."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from pathlib import Path


def parse_pure_or(gpr: str) -> list[str] | None:
    text = str(gpr or "").strip()
    if not text or re.search(r"\band\b|&", text, flags=re.I):
        return None
    # A single outer parenthesis around a pure OR expression is common in CarveMe.
    while text.startswith("(") and text.endswith(")"):
        depth = 0
        wraps_all = True
        for i, char in enumerate(text):
            if char == "(":
                depth += 1
            elif char == ")":
                depth -= 1
                if depth == 0 and i != len(text) - 1:
                    wraps_all = False
                    break
        if not wraps_all:
            break
        text = text[1:-1].strip()
    if "(" in text or ")" in text:
        return None
    genes = [token.strip() for token in re.split(r"\bor\b|\|", text, flags=re.I)]
    genes = [gene for gene in genes if gene]
    if len(genes) < 2:
        return None
    return list(dict.fromkeys(genes))


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_dump", required=True)
    parser.add_argument("--kofam", required=True)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--gene_pos_table", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--min_or", type=int, default=2)
    args = parser.parse_args()

    helper_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(helper_dir))
    from gene_mapping import (
        load_gff_records,
        load_gene_pos_table,
        load_kofam_index,
        resolve_gene_ids,
    )
    from gpr_evidence import decide_reaction, gene_ec_map

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    model_rows = read_tsv(Path(args.model_dump))
    reactions = []
    all_genes = []
    seen_genes = set()
    for row in model_rows:
        genes = parse_pure_or(row.get("gpr", ""))
        if genes is None or len(genes) < args.min_or:
            continue
        reactions.append((row, genes))
        for gene in genes:
            if gene not in seen_genes:
                seen_genes.add(gene)
                all_genes.append(gene)

    kofam = load_kofam_index(args.kofam)
    annotation_records = load_gff_records(args.gff)
    gene_pos = load_gene_pos_table(args.gene_pos_table)
    mapping, methods, unresolved = resolve_gene_ids(
        all_genes, kofam, annotation_records, gene_pos
    )
    gene_kos = {
        gene: list(kofam["by_id"].get(mapping.get(gene, ""), []))
        for gene in all_genes
    }
    gene_ecs = gene_ec_map({gene: kos for gene, kos in gene_kos.items() if kos})

    reaction_out = []
    gene_out = []
    category_counts = Counter()
    status_counts = Counter()
    for row, genes in reactions:
        rid = row["reaction_id"]
        reaction_ecs = {gene: gene_ecs.get(gene, set()) for gene in genes}
        decision = decide_reaction(
            rid,
            genes,
            reaction_ecs,
            gene_kos={gene: gene_kos.get(gene, []) for gene in genes},
            allow_add=False,
        )
        consensus = decision.consensus
        removals = decision.removals
        annotated = sum(bool(reaction_ecs.get(gene)) for gene in genes)
        no_ec = len(genes) - annotated
        if consensus.status == "confident" and removals:
            category = "confident_removal_candidate"
        elif consensus.status == "confident":
            category = "confident_no_change"
        elif consensus.status == "ambiguous":
            category = "ambiguous_needs_review"
        elif no_ec == len(genes):
            category = "no_ec_no_change"
        else:
            category = "insufficient_needs_review"
        category_counts[category] += 1
        status_counts[consensus.status] += 1
        clusters = ";".join(f"{ec}:{n}" for ec, n in consensus.clusters)
        reaction_out.append({
            "reaction_id": rid,
            "name": row.get("name", ""),
            "n_genes": len(genes),
            "n_mapped": sum(gene in mapping for gene in genes),
            "n_ec_annotated": annotated,
            "n_no_ec": no_ec,
            "status": consensus.status,
            "consensus_ec": consensus.ec or "",
            "confidence": f"{consensus.confidence:.4f}",
            "n_support": consensus.n_support,
            "clusters": clusters,
            "category": category,
            "n_rule_removals": len(removals),
            "rule_removals": ",".join(removals),
            "needs_llm": int(decision.needs_llm),
            "genes": ",".join(genes),
            "mapping_methods": ",".join(sorted({methods.get(g, "unmapped") for g in genes})),
            "gpr": row.get("gpr", ""),
        })
        for gene_decision in decision.decisions:
            gene = gene_decision.gene
            gene_out.append({
                "reaction_id": rid,
                "reaction_name": row.get("name", ""),
                "gene": gene,
                "gene_label": gene,
                "kofam_gene": mapping.get(gene, ""),
                "mapping_method": methods.get(gene, "unmapped"),
                "kos": ",".join(gene_decision.kos),
                "ecs": ",".join(gene_decision.ecs),
                "state": gene_decision.state,
                "action": gene_decision.action,
                "consensus_status": consensus.status,
                "consensus_ec": consensus.ec or "",
                "confidence": f"{consensus.confidence:.4f}",
                "n_support": consensus.n_support,
            })

    with (out_dir / "or2_reactions.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(reaction_out[0]) if reaction_out else ["reaction_id"])
        writer.writeheader()
        writer.writerows(reaction_out)
    with (out_dir / "or2_genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(gene_out[0]) if gene_out else ["reaction_id"])
        writer.writeheader()
        writer.writerows(gene_out)
    with (out_dir / "unresolved_genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene", "reason"])
        writer.writerows(sorted(unresolved.items()))

    summary = {
        "n_reactions_total": len(model_rows),
        "n_pure_or_ge_min": len(reaction_out),
        "n_model_genes": len(all_genes),
        "n_mapped_genes": len(mapping),
        "n_unresolved_genes": len(unresolved),
        "mapping_methods": dict(sorted(Counter(methods.values()).items())),
        "status_counts": dict(sorted(status_counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "n_reactions_with_rule_removals": sum(
            int(row["n_rule_removals"]) > 0 for row in reaction_out
        ),
        "n_rule_removed_gene_instances": sum(
            int(row["n_rule_removals"]) for row in reaction_out
        ),
        "n_reactions_needing_review": sum(
            row["category"].endswith("needs_review") for row in reaction_out
        ),
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
