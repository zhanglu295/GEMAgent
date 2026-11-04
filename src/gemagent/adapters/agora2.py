#!/usr/bin/env python3
"""Normalise an AGORA2 annotated SBML into a Cobra-readable, gapseq-style model.

- repairs invalid UTF-8 / control bytes that make the published XML malformed
- renames species  M_<cpd>__91__<comp>__93__  ->  M_<cpd>_<comp>0
- renames exchange R_EX_<cpd>__40__e__41__   ->  EX_<cpd>_e0

The mapping is a pure identifier rewrite; stoichiometry, bounds and the
objective are untouched.
"""
import re
import sys
from pathlib import Path


def normalize(src, dst):
    text = Path(src).read_bytes().decode("utf-8", errors="replace")
    text = text.replace("\ufffd", "?")
    text = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f]", "", text)

    counts = {"species": 0, "exchange": 0}

    def species(match):
        counts["species"] += 1
        return "M_%s_%s0" % (match.group(1), match.group(2))

    def exchange(match):
        counts["exchange"] += 1
        return "EX_%s_e0" % match.group(1)

    text = re.sub(r"\bM_([A-Za-z0-9_]+)__91__([a-z])__93__", species, text)
    text = re.sub(r"\bR_EX_([A-Za-z0-9_]+)__40__e__41__", exchange, text)
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    Path(dst).write_text(text, encoding="utf-8")
    return counts


if __name__ == "__main__":
    print(normalize(sys.argv[1], sys.argv[2]))
