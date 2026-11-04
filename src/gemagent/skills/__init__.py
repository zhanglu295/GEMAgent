"""Built-in review-skill registry (ordinary Python plugins, not Codex skills)."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SkillSpec:
    name: str
    description: str
    requires: tuple[str, ...]


SKILLS = {
    "biomass": SkillSpec(
        "biomass",
        "Audit biomass precursor requirements after gapseq draft.",
        ("draft", "medium", "kofam"),
    ),
    "gapfill-review": SkillSpec(
        "gapfill-review",
        "Review reactions added by official gapseq fill; exclude high-confidence bypasses and re-fill.",
        ("draft", "gapfilled", "medium", "exclude-list"),
    ),
    "gpr": SkillSpec(
        "gpr",
        "Review pure-OR branches with OR2-LLM-v1 and apply REMOVE decisions only.",
        ("gapfilled", "kofam"),
    ),
    "us-review-expand": SkillSpec(
        "us-review-expand",
        "Apply positive literature claims to add missing uptake/secretion capabilities.",
        ("gpr-checked", "claims", "exchange-map"),
    ),
    "us-review-prune": SkillSpec(
        "us-review-prune",
        "Audit existing positive calls with GEMAgent model evidence and organism context; KEEP/BLOCK/UNCERTAIN.",
        ("expanded-model", "draft", "gapfilled", "gprchecked", "claims", "context"),
    ),
    # compatibility alias: old name must not silently perform BLOCK decisions
    "us-review": SkillSpec(
        "us-review",
        "Compatibility alias for us-review-expand (positive literature additions only).",
        ("gpr-checked", "claims", "exchange-map"),
    ),
}


def get_skill(name: str) -> SkillSpec:
    if name not in SKILLS:
        raise KeyError(f"unknown skill: {name}")
    return SKILLS[name]
