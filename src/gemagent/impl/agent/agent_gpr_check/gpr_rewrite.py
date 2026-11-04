"""Pure-OR GPR rewriting used by both the cobra and SBML-text appliers."""
from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple


def strip_outer_parens(text: str) -> str:
    """Remove parentheses that wrap the entire expression (one or more pairs)."""
    s = str(text).strip()
    while s.startswith("(") and s.endswith(")"):
        depth = 0
        wraps = True
        for i, ch in enumerate(s):
            if ch == "(":
                depth += 1
            elif ch == ")":
                depth -= 1
                if depth == 0 and i < len(s) - 1:
                    wraps = False
                    break
        if not wraps:
            break
        s = s[1:-1].strip()
    return s


def gene_alias(gene: str) -> str:
    g = str(gene).strip()
    return g[2:] if g.startswith("G_") else g


def parse_isozyme(gpr: str) -> Optional[List[str]]:
    """Return the ordered gene list for a pure OR rule, else None."""
    s = strip_outer_parens(str(gpr or "").strip())
    if not s:
        return None
    if re.search(r"\band\b", s, re.I):
        return None
    genes = []
    for part in re.split(r"\s+or\s+", s, flags=re.I):
        token = strip_outer_parens(part.strip())
        if not token or "(" in token or ")" in token:
            return None
        genes.append(token)
    return genes if genes else None


def action_for(gene: str, actions: Dict[str, str]) -> str:
    """Look up an action tolerating the G_ prefix used by CarveMe."""
    base = gene_alias(gene)
    for key in (str(gene).strip(), base, "G_" + base):
        if key in actions:
            return (actions[key] or "").upper()
    lower = {gene_alias(k).lower(): v for k, v in actions.items()}
    return (lower.get(base.lower()) or "").upper()


def apply_actions(gpr: str, actions: Dict[str, str]) -> Tuple[str, List[str], List[str], bool]:
    """Apply KEEP/REMOVE/ADD to a pure OR rule.

    Returns (new_gpr, removed, added, skipped_all_removed).
    If every gene would be removed, the original rule is kept.
    """
    genes = parse_isozyme(gpr)
    if not genes:
        return gpr, [], [], False

    remove = {g for g in genes if action_for(g, actions) == "REMOVE"}
    keep = [g for g in genes if g not in remove]
    removed = [g for g in genes if g in remove]

    if not keep:
        return gpr, [], [], True

    use_g_prefix = all(g.startswith("G_") for g in genes)
    add = []
    for gene, action in actions.items():
        if (action or "").upper() != "ADD":
            continue
        base = gene_alias(gene)
        if any(gene_alias(g) == base for g in genes):
            continue
        add.append("G_" + base if use_g_prefix else base)

    new_gpr = " or ".join(keep + add)
    return new_gpr, removed, add, False
