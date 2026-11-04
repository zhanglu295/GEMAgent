#!/usr/bin/env python3
"""EC-evidence layer for GPR review.

The review used to be a deletion-only isozyme pruner driven by an LLM prior.
This module replaces that prior with a reaction-level EC consensus computed
from the KOfam annotations of the genes already attached to the reaction.

The rules are deliberately conservative:

* KEEP is the default action.
* REMOVE requires a full (4-level) EC conflict against a confident consensus.
* ADD requires a confident consensus, a candidate gene with an exact EC match,
  no current gene supporting that EC, and at most `max_add` additions.
* When the evidence is insufficient or ambiguous, no edit is made and the
  reaction is flagged for optional LLM adjudication.

No curated model is used anywhere in this module.
"""

from __future__ import annotations

import re
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple

EC_TOKEN = re.compile(r"\d+\.\d+\.\d+\.(?:\d+|-|n)")
EC_BLOCK = re.compile(r"\[EC:([^\]]+)\]")


def parse_ecs_from_definition(definition: str) -> Set[str]:
    """Extract EC numbers from a KOfam KO definition."""
    ecs: Set[str] = set()
    for block in EC_BLOCK.findall(definition or ""):
        for token in re.split(r"[\s,;]+", block):
            token = token.strip()
            if EC_TOKEN.fullmatch(token):
                ecs.add(token)
    return ecs


def is_full_ec(ec: str) -> bool:
    """True for a complete four-level EC number (no '-' or 'n')."""
    if not ec:
        return False
    parts = ec.split(".")
    return len(parts) == 4 and all(p.isdigit() for p in parts)


def ec_prefix(ec: str, levels: int = 3) -> str:
    return ".".join(ec.split(".")[:levels]) + "."


def gene_ec_map(gene_kos: Dict[str, Sequence[Tuple[str, str]]]) -> Dict[str, Set[str]]:
    """gene -> union of EC numbers across its best-hit KOs."""
    out: Dict[str, Set[str]] = {}
    for gene, kos in gene_kos.items():
        ecs: Set[str] = set()
        for _ko, definition in kos:
            ecs |= parse_ecs_from_definition(definition)
        out[str(gene)] = ecs
    return out


def gene_weight(ecs: Iterable[str], max_ec_per_gene: int = 3) -> float:
    """Weight a gene by annotation specificity; promiscuous genes do not vote."""
    ecs = list(ecs)
    n = len(ecs)
    if n == 0 or n > max_ec_per_gene:
        return 0.0
    return 1.0 / n


@dataclass
class Consensus:
    status: str  # confident | ambiguous | insufficient
    ec: Optional[str]
    confidence: float
    n_support: int
    support: Dict[str, float] = field(default_factory=dict)
    clusters: List[Tuple[str, int]] = field(default_factory=list)

    @property
    def confident(self) -> bool:
        return self.status == "confident"


def ec_consensus(
    gene_ecs: Dict[str, Set[str]],
    min_support: int = 2,
    min_support_frac: float = 0.60,
    max_ec_per_gene: int = 3,
) -> Consensus:
    """Weighted vote over the full EC numbers carried by the current genes."""
    weights: Dict[str, float] = {}
    total = 0.0
    for gene, ecs in gene_ecs.items():
        w = gene_weight(ecs, max_ec_per_gene)
        weights[gene] = w
        total += w

    support: Dict[str, float] = defaultdict(float)
    support_n: Dict[str, int] = defaultdict(int)
    for gene, ecs in gene_ecs.items():
        w = weights[gene]
        if w <= 0:
            continue
        for ec in ecs:
            if is_full_ec(ec):
                support[ec] += w
                support_n[ec] += 1

    if not support or total <= 0:
        return Consensus(status="insufficient", ec=None, confidence=0.0, n_support=0)

    ranked = sorted(support.items(), key=lambda kv: (kv[1], support_n[kv[0]]), reverse=True)
    best_ec, best_weight = ranked[0]
    confidence = best_weight / total
    clusters = [(ec, support_n[ec]) for ec, _w in ranked if support_n[ec] > 0]

    if support_n[best_ec] >= min_support and confidence >= min_support_frac:
        status = "confident"
    elif len([c for c in clusters if c[1] >= 1]) >= 2:
        status = "ambiguous"
    else:
        status = "insufficient"

    return Consensus(
        status=status,
        ec=best_ec,
        confidence=confidence,
        n_support=support_n[best_ec],
        support=dict(support),
        clusters=clusters,
    )


def classify_gene(ecs: Set[str], consensus: Consensus) -> str:
    """support | partial | conflict | unknown relative to the consensus EC."""
    if not consensus.confident or not consensus.ec:
        return "unknown"
    if not ecs:
        return "unknown"
    if consensus.ec in ecs:
        return "support"
    prefix = ec_prefix(consensus.ec)
    if any(e.startswith(prefix) for e in ecs):
        return "partial"
    return "conflict"


@dataclass
class GeneDecision:
    gene: str
    action: str  # KEEP | REMOVE | ADD
    state: str
    ecs: List[str]
    kos: List[str]
    evidence_ec: str = ""
    evidence_ko: str = ""
    source: str = "rule"


@dataclass
class ReactionDecision:
    reaction_id: str
    genes: List[str]
    consensus: Consensus
    decisions: List[GeneDecision]
    needs_llm: bool
    notes: str = ""

    @property
    def removals(self) -> List[str]:
        return [d.gene for d in self.decisions if d.action == "REMOVE"]

    @property
    def additions(self) -> List[str]:
        return [d.gene for d in self.decisions if d.action == "ADD"]


def decide_reaction(
    reaction_id: str,
    genes: Sequence[str],
    gene_ecs: Dict[str, Set[str]],
    gene_kos: Dict[str, Sequence[Tuple[str, str]]] | None = None,
    candidate_ecs: Dict[str, Set[str]] | None = None,
    candidate_kos: Dict[str, Sequence[Tuple[str, str]]] | None = None,
    allow_add: bool = True,
    max_add: int = 2,
    min_support: int = 2,
    min_support_frac: float = 0.60,
    max_ec_per_gene: int = 3,
) -> ReactionDecision:
    """Apply the conservative evidence rules to one pure-OR reaction."""
    genes = [str(g) for g in genes]
    gene_kos = gene_kos or {}
    gene_ecs = {str(g): set(gene_ecs.get(g, set())) for g in genes}
    consensus = ec_consensus(gene_ecs, min_support, min_support_frac, max_ec_per_gene)

    decisions: List[GeneDecision] = []
    states = {g: classify_gene(gene_ecs[g], consensus) for g in genes}

    # Default: keep everything.
    for g in genes:
        kos = [ko for ko, _d in gene_kos.get(g, [])]
        decisions.append(
            GeneDecision(
                gene=g,
                action="KEEP",
                state=states[g],
                ecs=sorted(gene_ecs[g]),
                kos=kos,
            )
        )

    if consensus.confident:
        conflicting = [g for g in genes if states[g] == "conflict"]
        # Never empty a GPR: require at least one non-conflicting gene.
        survivors = [g for g in genes if g not in conflicting]
        if conflicting and survivors:
            by_gene = {d.gene: d for d in decisions}
            for g in conflicting:
                by_gene[g].action = "REMOVE"
                conflicting_ec = next(
                    (e for e in sorted(gene_ecs[g]) if is_full_ec(e)), ""
                )
                by_gene[g].evidence_ec = conflicting_ec
                by_gene[g].evidence_ko = by_gene[g].kos[0] if by_gene[g].kos else ""

        # Add only when nothing in the current GPR supports the consensus EC.
        if allow_add and max_add > 0 and candidate_ecs and not any(
            states[g] == "support" for g in genes
        ):
            added = 0
            for cand in sorted(candidate_ecs):
                if cand in genes or added >= max_add:
                    continue
                cecs = set(candidate_ecs[cand])
                if consensus.ec not in cecs or len(cecs) > 2:
                    continue
                kos = [ko for ko, _d in (candidate_kos or {}).get(cand, [])]
                decisions.append(
                    GeneDecision(
                        gene=cand,
                        action="ADD",
                        state="support",
                        ecs=sorted(cecs),
                        kos=kos,
                        evidence_ec=consensus.ec or "",
                        evidence_ko=kos[0] if kos else "",
                    )
                )
                added += 1

    has_ec = [bool(gene_ecs[g]) for g in genes]
    mixed_annotation = any(has_ec) and not all(has_ec)
    needs_llm = consensus.status == "ambiguous" or (
        consensus.status == "insufficient" and mixed_annotation
    )
    notes = ""
    if consensus.status == "insufficient":
        notes = ("mixed EC/no-EC annotation; adjudication required" if mixed_annotation
                 else "no EC evidence; no edit")
    elif consensus.status == "ambiguous":
        notes = "multiple EC clusters; adjudication required"

    return ReactionDecision(
        reaction_id=reaction_id,
        genes=genes,
        consensus=consensus,
        decisions=decisions,
        needs_llm=needs_llm,
        notes=notes,
    )
