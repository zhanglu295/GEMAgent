#!/usr/bin/env python3
"""Single source of truth for every filesystem / external-binary location.

Before this module existed each component re-derived the project root with its
own ``dirname()`` chain and hardcoded absolute installation paths for Python,
Rscript, and gapseq. That made the tree non-relocatable and let
two different SEED databases (``code/dat`` = the installed gapseq's own copy,
``database/`` = GEMAgent's curated copy) drift apart unnoticed -- they differ in
the name of their status column, ``gapseq.status`` vs ``gemagent.status``,
which used to be hardcoded at the call sites.

Everything resolves relative to this file, with environment overrides first:

    GEMAGENT_PYTHON / PYTHON       python interpreter   (default: sys.executable)
    GEMAGENT_RSCRIPT / RSCRIPT     Rscript binary       (default: sibling of python)
    GAPSEQ_BIN                     gapseq CLI           (default: which, then conda envs)

Database roles (they are NOT interchangeable):

    code/dat/      the *installed* gapseq's native ``dat/``, reached through
                   the symlink ``code/dat -> <gapseq env>/share/gapseq/dat``
                   (currently bioconda gapseq 1.4.0 -- check with
                   ``gapseq -v`` / ``readlink code/dat``). The gapfilling R
                   scripts (gf.suite.R and friends) read this one through
                   ``file.path(script.dir, "..", "dat", ...)``.
    database/      GEMAgent's curated copy. Partly gapseq-2.1-derived, and it
                   also holds files gapseq's dat does not ship
                   (mimedb/mimedb_microbes_v2.csv,
                   media/medium_prediction_rules.tsv, metabolite/diffusion_mets.tsv).
                   The biomass / medium / prompt code reads this one, through
                   either ``paths.db(...)`` or an explicit ``database/`` path in
                   R (replace_biomass.R, predict_medium.R).

                   ``database/taxonomy/`` is the offline naming cache (GTDB
                   r214 lineage + NCBI organism names); it is read only through
                   :mod:`taxonomy`, which also documents why the GTDB and NCBI
                   names must not be interchanged.

``database/`` deliberately duplicates a few files that also exist in ``code/dat``
(``biomass/*.json``, ``media/ALLmed.csv``, ``metabolite/seed_metabolites_edited.tsv``)
because those two consumers read them from different roots. Everything else that
used to sit here unused (``taxonomy.tsv``, ``media/MM_glu.csv``,
``biomass/excluded_reactions.tsv``, ``reaction/mnxref_*``, ``transporter/subex.tsv``)
was moved to ``_unused_20260919/`` -- the gapfilling layer reads its copies from
``code/dat``.

The two do overlap on ``seed_reactions_corrected.tsv``. Use
:func:`seed_reactions_tsv` to obtain it together with its actual status-column
name, and :func:`status_column` to read that column defensively.

DO NOT point the gapfilling layer at ``database/reaction/seed_reactions_corrected.tsv``.
Tried on 2026-09-20 and reverted: that table marks 7 reactions as ``removed``
that ``code/dat`` has as approved/corrected -- including BOTH NAD(P)
transhydrogenases (rxn10125, rxn27378) and three ubiquinol oxidases. Without
NADH<->NADPH interconversion the Gram-positive biomass template, which requires
both redox currencies, becomes infeasible; gapfilling then aborts with "cannot
produce enough target even when all candidate reactions are added" and the
model keeps growth = 0 (observed on cper GCF_000013285.1).
"""
from __future__ import annotations

import os
import shutil
import sys

# ---------------------------------------------------------------------------
# Tree layout
# ---------------------------------------------------------------------------
CODE_DIR = os.path.dirname(os.path.abspath(__file__))
PACKAGE_DIR = os.path.dirname(CODE_DIR)
SOURCE_ROOT = os.path.dirname(PACKAGE_DIR)
PROJECT_ROOT = os.path.dirname(SOURCE_ROOT)

#: gapseq 2.1 native dat (symlink). Read by the gapfilling R scripts.
DAT_DIR = os.environ.get("GEMAGENT_GAPSEQ_DAT_DIR", os.path.join(CODE_DIR, "dat"))
#: GEMAgent's curated database. Read by the biomass / medium / prompt code.
DB_DIR = os.environ.get("GEMAGENT_DB_DIR", os.path.join(PROJECT_ROOT, "database"))


def db(*parts: str) -> str:
    """Path inside GEMAgent's ``database/`` (e.g. ``db("biomass", "x.json")``)."""
    return os.path.join(DB_DIR, *parts)


def dat(*parts: str) -> str:
    """Path inside gapseq 2.1's ``dat/`` (e.g. ``dat("seed_reactions_corrected.tsv")``)."""
    return os.path.join(DAT_DIR, *parts)


# ---------------------------------------------------------------------------
# External binaries
# ---------------------------------------------------------------------------
def _env(*names: str) -> str:
    for n in names:
        v = os.environ.get(n, "").strip()
        if v:
            return v
    return ""


def _conda_env_bin(name: str) -> str:
    """Look for ``name`` in any conda env under this interpreter's conda root.

    The cluster layout keeps gapseq in its own env (``envs/gapseq/bin``), so a
    plain ``shutil.which`` from inside the gemagent env will not find it. Search
    sibling envs instead of hardcoding an absolute path.
    """
    env_bin = os.path.dirname(os.path.abspath(sys.executable))      # .../envs/<env>/bin
    env_dir = os.path.dirname(env_bin)                             # .../envs/<env>
    envs_dir = os.path.dirname(env_dir)                            # .../envs
    if not os.path.isdir(envs_dir):
        return ""
    for entry in sorted(os.listdir(envs_dir)):
        cand = os.path.join(envs_dir, entry, "bin", name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    return ""


def resolve_bin(name: str, *env_names: str, sibling_ok: bool = False) -> str:
    """Resolve an external executable without hardcoding an absolute path.

    Order: environment override -> a sibling of the running interpreter
    (``sibling_ok``; python and Rscript share a conda env) -> ``PATH`` ->
    any conda env under the conda root.
    """
    got = _env(*env_names)
    if got:
        return got
    if sibling_ok:
        cand = os.path.join(os.path.dirname(os.path.abspath(sys.executable)), name)
        if os.path.isfile(cand) and os.access(cand, os.X_OK):
            return cand
    found = shutil.which(name)
    if found:
        return found
    return _conda_env_bin(name)


#: Python interpreter running the pipeline.
PYTHON = _env("GEMAGENT_PYTHON", "PYTHON") or os.path.abspath(sys.executable)
#: Rscript from the same conda env as the interpreter, else PATH, else any env.
RSCRIPT = resolve_bin("Rscript", "GEMAGENT_RSCRIPT", "RSCRIPT", sibling_ok=True)
#: Native gapseq CLI (its own env).
GAPSEQ_BIN = resolve_bin("gapseq", "GAPSEQ_BIN")


# --- KOFamScan toolchain (gene functional annotation, step_annotate_kofam) ----
# `prodigal` and `exec_annotation` live in *different* conda envs (gapseq and
# deepmm_eval) and are therefore looked up across envs, never hardcoded. The HMM
# profile and ko_list are data files that ship with neither, so they are taken
# from KOFAM_DIR (default ``~/kofam_local``), which is the one path an operator
# may need to move.
KOFAM_DIR = _env("GEMAGENT_KOFAM_DIR", "KOFAM_DIR") or os.path.join(
    os.path.expanduser("~"), "kofam_local")
PRODIGAL = resolve_bin("prodigal", "GEMAGENT_PRODIGAL", "PRODIGAL")
EXEC_ANNOTATION = resolve_bin("exec_annotation", "GEMAGENT_EXEC_ANNOTATION",
                             "EXEC_ANNOTATION")
KOFAM_PROFILE = _env("KOFAM_PROFILE") or os.path.join(
    KOFAM_DIR, "profiles", "prokaryote.hal")
KOFAM_KO_LIST = _env("KOFAM_KO_LIST") or os.path.join(KOFAM_DIR, "ko_list")

#: Env var -> resolved value, for handing the whole toolchain to annotate_kofam.sh
#: (which reads exactly these names, and falls back to its own defaults when run
#: standalone).
KOFAM_ENV = {
    "PRODIGAL": PRODIGAL,
    "EXEC_ANNOTATION": EXEC_ANNOTATION,
    "KOFAM_PROFILE": KOFAM_PROFILE,
    "KOFAM_KO_LIST": KOFAM_KO_LIST,
}


def missing_kofam_tools() -> list[str]:
    """Human-readable list of the annotation inputs that are absent."""
    out = []
    for label, path in (("PRODIGAL", PRODIGAL), ("EXEC_ANNOTATION", EXEC_ANNOTATION),
                        ("KOFAM_PROFILE", KOFAM_PROFILE),
                        ("KOFAM_KO_LIST", KOFAM_KO_LIST)):
        if not path or not os.path.exists(path):
            out.append(f"{label}={path or '<unresolved>'}")
    return out


# ---------------------------------------------------------------------------
# The two SEED reaction tables
# ---------------------------------------------------------------------------
#: Column name that marks a reaction as suppressed, per database flavour.
STATUS_COLUMNS = ("gapseq.status", "gemagent.status", "status")


def status_column(header) -> str:
    """Return whichever suppression column *header* actually carries.

    ``code/dat`` ships ``gapseq.status`` (gapseq 2.1) while ``database/`` ships
    ``gemagent.status`` (GEMAgent's curated copy); the two trees are otherwise
    column-for-column identical. Hardcoding either name makes the reader break
    whenever it is pointed at the other tree.
    """
    for col in STATUS_COLUMNS:
        if col in header:
            return col
    return ""


def seed_reactions_tsv(prefer: str = "dat") -> str:
    """Locate ``seed_reactions_corrected.tsv``; argsort-free with a clear order.

    ``prefer="dat"`` (default) picks gapseq's own shipped copy, which is the
    table the gapfilling layer actually reads (``code/gapfilling/*.R`` resolve
    ``../dat/``), so any status decision taken from it stays consistent with
    what gapfilling did. ``prefer="database"`` falls back to GEMAgent's curated
    copy.

    Do NOT flip this default back to ``"database"``: see the long note in the
    module docstring -- the curated copy has both NAD(P) transhydrogenases
    marked ``removed``, which makes every Gram-positive biomass template
    infeasible and silently yields zero-growth models.
    """
    order = (dat("seed_reactions_corrected.tsv"),
             os.path.join(DB_DIR, "reaction", "seed_reactions_corrected.tsv"))
    if prefer == "database":
        order = order[::-1]
    for p in order:
        if os.path.exists(p):
            return p
    return order[0]


def seed_reactions_path() -> tuple[str, str]:
    """Return ``(path, status_column)`` for the SEED reaction table.

    The column is detected from the real header, so callers never have to know
    which of the two database flavours they were handed.
    """
    path = seed_reactions_tsv()
    if not os.path.exists(path):
        return path, ""
    with open(path, encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
    return path, status_column(header)


def taxonomy_table(*parts: str) -> str:
    """Path inside ``database/taxonomy/`` (the offline naming cache).

    Only :mod:`taxonomy` should read these -- it owns the fallback chain and
    the GTDB-vs-NCBI distinction. Kept here so the four table locations are
    declared in one place alongside the rest of the data layout.
    """
    return db("taxonomy", *parts)


def describe() -> str:
    """One-line-per-item dump, handy for --verbose / bug reports."""
    lines = [
        f"PROJECT_ROOT {PROJECT_ROOT}",
        f"CODE_DIR     {CODE_DIR}",
        f"DAT_DIR      {DAT_DIR}{'' if os.path.exists(DAT_DIR) else '  (MISSING)'}",
        f"DB_DIR       {DB_DIR}{'' if os.path.exists(DB_DIR) else '  (MISSING)'}",
        f"PYTHON       {PYTHON}",
        f"RSCRIPT      {RSCRIPT or '<NOT FOUND>'}",
        f"GAPSEQ_BIN   {GAPSEQ_BIN or '<NOT FOUND>'}",
        f"seed_rxns    {seed_reactions_tsv()}",
        f"taxonomy/    {taxonomy_table()}"
        f"{'' if os.path.isdir(taxonomy_table()) else '  (MISSING)'}",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    print(describe())
