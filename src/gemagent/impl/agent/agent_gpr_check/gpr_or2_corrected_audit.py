#!/usr/bin/env python3
"""Audit all pure-OR GPRs with correct KOfam gene-ID resolution.

This script intentionally avoids COBRA so it can run in the lightweight system
Python environment. It parses SBML directly and reuses the project's existing
gene_mapping.py and gpr_evidence.py policies.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from pathlib import Path


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def children(node):
    return list(node) if node is not None else []


def first_child(node):
    kids = children(node)
    return kids[0] if kids else None


def child_by_local(node, name):
    for child in children(node):
        if local_name(child.tag) == name:
            return child
    return None


def attr_by_local(node, name):
    for key, value in node.attrib.items():
        if local_name(key) == name:
            return value
    return None


def descendants_by_local(node, name):
    return [x for x in node.iter() if local_name(x.tag) == name]


def parse_pure_or(association):
    """Return ordered geneProduct ids for a pure OR association, else None."""
    if association is None:
        return None
    root = first_child(association)
    if root is None or local_name(root.tag) != "or":
        return None
    genes = []
    for child in children(root):
        if local_name(child.tag) != "geneProductRef":
            return None
        gene = attr_by_local(child, "geneProduct") or attr_by_local(child, "geneProductId")
        if gene:
            genes.append(gene)
    return genes


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model", required=True)
    parser.add_argument("--kofam", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--project_root", default=".")
    parser.add_argument("--min_or", type=int, default=2)
    parser.add_argument("--gff", default="")
    args = parser.parse_args()

    helper_dir = Path(__file__).resolve().parent
    sys.path.insert(0, str(helper_dir))

    from gene_mapping import load_kofam_index, load_gff_records, resolve_gene_ids
    from gpr_evidence import decide_reaction, gene_ec_map

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    kofam = load_kofam_index(args.kofam)
    root = ET.parse(args.model).getroot()

    gene_labels = {}
    for node in descendants_by_local(root, "geneProduct"):
        gid = node.attrib.get("id") or attr_by_local(node, "id")
        if gid:
            gene_labels[gid] = attr_by_local(node, "label") or attr_by_local(node, "name") or ""

    reactions = descendants_by_local(root, "reaction")
    model_gene_ids = sorted(gene_labels)
    annotation_records = load_gff_records(args.gff) if args.gff else []
    mapping, methods, unresolved = resolve_gene_ids(
        model_gene_ids, kofam, annotation_records=annotation_records, gene_pos={}
    )
    gene_kos = {
        gid: list(kofam["by_id"].get(mapping.get(gid, ""), []))
        for gid in model_gene_ids
    }
    gene_ecs = gene_ec_map({gid: kos for gid, kos in gene_kos.items() if kos})

    reaction_rows = []
    gene_rows = []
    skipped_non_or = 0
    skipped_lt_min = 0

    for reaction in reactions:
        rid = reaction.attrib.get("id", "")
        name_node = child_by_local(reaction, "name")
        name = name_node.text.strip() if name_node is not None and name_node.text else ""
        assoc = child_by_local(reaction, "geneProductAssociation")
        genes = parse_pure_or(assoc)
        if genes is None:
            skipped_non_or += 1
            continue
        if len(genes) < args.min_or:
            skipped_lt_min += 1
            continue

        reaction_gene_kos = {g: gene_kos.get(g, []) for g in genes}
        reaction_gene_ecs = {g: gene_ecs.get(g, set()) for g in genes}
        decision = decide_reaction(
            rid,
            genes,
            reaction_gene_ecs,
            gene_kos=reaction_gene_kos,
            allow_add=False,
        )
        consensus = decision.consensus
        removals = decision.removals
        mapped = sum(1 for g in genes if g in mapping)
        annotated = sum(1 for g in genes if reaction_gene_ecs.get(g))
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

        clusters = ";".join(f"{ec}:{n}" for ec, n in consensus.clusters)
        reaction_rows.append({
            "reaction_id": rid,
            "name": name,
            "n_genes": len(genes),
            "n_mapped": mapped,
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
        })
        for gene_decision in decision.decisions:
            gene = gene_decision.gene
            gene_rows.append({
                "reaction_id": rid,
                "reaction_name": name,
                "gene": gene,
                "gene_label": gene_labels.get(gene, ""),
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
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(reaction_rows[0]) if reaction_rows else ["reaction_id"])
        writer.writeheader()
        writer.writerows(reaction_rows)

    with (out_dir / "or2_genes.tsv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=list(gene_rows[0]) if gene_rows else ["reaction_id"])
        writer.writeheader()
        writer.writerows(gene_rows)

    category_counts = Counter(row["category"] for row in reaction_rows)
    status_counts = Counter(row["status"] for row in reaction_rows)
    method_counts = Counter(methods.values())
    summary = {
        "model": os.path.abspath(args.model),
        "kofam": os.path.abspath(args.kofam),
        "n_reactions_total": len(reactions),
        "n_pure_or_ge_min": len(reaction_rows),
        "n_skipped_non_pure_or": skipped_non_or,
        "n_skipped_below_min_or": skipped_lt_min,
        "n_model_genes": len(model_gene_ids),
        "n_mapped_genes": len(mapping),
        "n_unresolved_genes": len(unresolved),
        "mapping_methods": dict(sorted(method_counts.items())),
        "status_counts": dict(sorted(status_counts.items())),
        "category_counts": dict(sorted(category_counts.items())),
        "n_reactions_with_rule_removals": sum(
            1 for row in reaction_rows if int(row["n_rule_removals"]) > 0
        ),
        "n_rule_removed_gene_instances": sum(
            int(row["n_rule_removals"]) for row in reaction_rows
        ),
        "n_reactions_needing_review": sum(
            1 for row in reaction_rows if row["category"].endswith("needs_review")
        ),
    }
    (out_dir / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    unresolved_path = out_dir / "unresolved_genes.tsv"
    with unresolved_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["gene", "reason"])
        for gene, reason in sorted(unresolved.items()):
            writer.writerow([gene, reason])

    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
