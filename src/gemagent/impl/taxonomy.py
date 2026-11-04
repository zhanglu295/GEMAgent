#!/usr/bin/env python3
"""Single source of truth for organism *naming* and *taxonomy*.

Two fields, two audiences -- they are NOT interchangeable:

    gtdb_lineage   GTDB r214 lineage, e.g.
                   ``d__Bacteria;p__Bacillota;c__Bacilli;o__Bacillales;
                   f__Bacillaceae;g__Bacillus;s__Bacillus subtilis``
                   -> ``taxonomy_tree`` prompt context for
                   biomass_check / review_gapfill_added / gpr_check.

    org_name       the **NCBI** organism name, e.g.
                   ``Clostridium perfringens ATCC 13124``
                   -> org_prompt_gen.py -> MIMEdb profile, and literature
                   queries in the deep-research module.

Why they must differ
--------------------
GTDB renames species that MIMEdb and PubMed do not follow. Measured against
``database/mimedb/mimedb_microbes_v2.csv`` (2648 rows; ``get_mimedb_info``
matches with ``str.contains`` on the ``species`` column):

    NCBI name                  rows   GTDB name                  rows
    Clostridium perfringens       4   Sarcina perfringens           0
    Eubacterium rectale           3   Agathobacter rectalis         0
    Clostridium acetobutylicum    2   Clostridium_S acetobutylicum  0
    Methanobrevibacter smithii    2   Methanobrevibacter_A smithii  0

Feeding the GTDB name to ``org_prompt_gen`` therefore silently drops the
physiology context of every review ("No database records found in MimeDB"),
and the same string retrieves nothing from literature search. GTDB also
decorates genera with suffixes (``_S``/``_A``/``_I``) that exist in no
database. So: GTDB drives ``taxonomy_tree``, NCBI drives ``org_name``.

Resolution order (first hit wins, recorded in the ``source`` column)
-------------------------------------------------------------------
    gtdb_lineage
        1. GTDB r214 bulk table      database/taxonomy/gtdb_all_r214.tsv
           (or the two ``*_taxonomy_r214.tsv.gz`` it is built from)
        2. GTDB-TK run               only with --gtdb_fallback (slow: pplacer)
    org_name
        1. NCBI RefSeq bulk table    database/taxonomy/ncbi_accession2name.tsv
           (or assembly_summary_refseq.txt it is built from)
        2. NCBI Datasets API         online, ~1 s/accession, cached on disk
        3. GTDB species, suffix-stripped as a last resort

Why the bulk tables and not a GTDB-TK run per genome
----------------------------------------------------
``gtdbtk`` ships only a **representative-genome** taxonomy
(``GTDBTK_DATA_PATH/taxonomy/gtdb_taxonomy.tsv`` = 85,205 rows). Measured on
the 26 FermentationTest genomes it resolves 16/26; E. coli K-12
(GCF_000005845.2) and 9 other named strains are absent. The 4.3 MB release
table (``bac120_taxonomy_r214.tsv.gz`` + ``ar53_taxonomy_r214.tsv.gz``,
402,709 rows) resolves **25/26** in milliseconds, versus 10-30 min/genome for
pplacer. Only ``GCF_000195895.1`` (Methanosarcina barkeri str. Fusaro) is
genuinely absent from r214.

Accession matching must be normalised: GTDB rewrites the accession prefix and
sometimes the assembly letter, e.g. the same genome is ``GCF_000020605.1`` in
RefSeq but ``GB_GCA_000020605.1`` in GTDB. Match on the numeric body.
"""
from __future__ import annotations

import gzip
import os
import re
import urllib.error
import urllib.request

import paths

# ---------------------------------------------------------------------------
# Tables and output contract
# ---------------------------------------------------------------------------
GTDB_MERGED = paths.taxonomy_table("gtdb_all_r214.tsv")
GTDB_PARTS = (paths.taxonomy_table("bac120_taxonomy_r214.tsv.gz"),
              paths.taxonomy_table("ar53_taxonomy_r214.tsv.gz"))
NCBI_SLIM = paths.taxonomy_table("ncbi_accession2name.tsv")
NCBI_BULK = paths.taxonomy_table("assembly_summary_refseq.txt")
API_CACHE = paths.taxonomy_table("ncbi_api_cache.tsv")

#: Written under ``<output_dir>/<OUT_SUBDIR>/<OUT_FILE>`` -- the only file the
#: review steps read, so a run's naming is reproducible and auditable.
OUT_SUBDIR = "taxonomy"
OUT_FILE = "taxonomy.tsv"
OUT_COLUMNS = ("genome_id", "accession", "gtdb_lineage", "gtdb_species",
               "org_name", "species_name", "source")

#: ``GCF_000013285.1`` / ``GCA_000013285`` -- RefSeq/GenBank assembly accession.
#: The trailing boundary is deliberately absent: a ``\b`` cannot hold between the
#: version digits and the following ``_`` (as in ``GCF_000013285.1_ASM1328v1``),
#: so the optional version group would backtrack away and silently drop ``.1`` --
#: and the NCBI table is keyed on the versioned accession. The lookbehind keeps
#: ``foo_GCF_...`` from being clipped mid-token.
ACC_RE = re.compile(r"(?<![A-Za-z0-9])(GC[AF]_\d{5,9}(?:\.\d+)?)")
_NUM_RE = re.compile(r"_(\d{5,9})(?:\.(\d+))?")

DATASETS_API = ("https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/"
                "{acc}/dataset_report")


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------
class TaxonomyResult:
    """Resolved naming for one genome."""

    __slots__ = ("genome_id", "accession", "gtdb_lineage", "gtdb_species",
                 "org_name", "species_name", "source")

    def __init__(self, genome_id="", accession="", gtdb_lineage="",
                 gtdb_species="", org_name="", species_name="", source=""):
        self.genome_id = genome_id
        self.accession = accession
        self.gtdb_lineage = gtdb_lineage
        self.gtdb_species = gtdb_species
        self.org_name = org_name
        self.species_name = species_name
        self.source = source

    def as_row(self) -> list[str]:
        return [getattr(self, c) for c in OUT_COLUMNS]

    def __repr__(self) -> str:
        return (f"TaxonomyResult(acc={self.accession!r}, "
                f"org_name={self.org_name!r}, source={self.source!r}, "
                f"gtdb_species={self.gtdb_species!r})")

    def describe(self) -> str:
        return "\n".join([
            f"accession    {self.accession or '<none>'}",
            f"org_name     {self.org_name or '<none>'}   (NCBI, -> literature)",
            f"species_name {self.species_name or '<none>'}   (-> MIMEdb lookup)",
            f"gtdb_species {self.gtdb_species or '<none>'}",
            f"gtdb_lineage {self.gtdb_lineage or '<none>'}",
            f"source       {self.source or '<none>'}",
        ])


# ---------------------------------------------------------------------------
# Accession handling
# ---------------------------------------------------------------------------
def normalise(acc: str) -> str:
    """Numeric body of an accession, e.g. ``GCF_000020605.1`` -> ``000020605``.

    GTDB stores the same genome under a different prefix/letter
    (``GB_GCA_000020605.1``), so the letters cannot be compared.
    """
    m = _NUM_RE.search(acc or "")
    return m.group(1) if m else ""


def extract_accession(*texts: str) -> str:
    """First ``GC[AF]_...`` accession found across *texts*, in order."""
    for t in texts:
        if not t:
            continue
        m = ACC_RE.search(t)
        if m:
            return m.group(1)
    return ""


def _fasta_first_header(genome_file: str) -> str:
    """``>`` line of a (possibly gzipped) FASTA, '' when unavailable."""
    if not genome_file or not os.path.exists(genome_file):
        return ""
    opener = gzip.open if genome_file.endswith(".gz") else open
    try:
        with opener(genome_file, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith(">"):
                    return line[1:].strip()
                if line.strip():
                    break
    except OSError:
        return ""
    return ""


# ---------------------------------------------------------------------------
# Lazy table readers (module-level cache: load once per process)
# ---------------------------------------------------------------------------
_gtdb_cache: dict[str, tuple[str, str]] | None = None
_ncbi_cache: dict[str, str] | None = None
_ncbi_body_cache: dict[str, str] | None = None


def _iter_gtdb_lines():
    if os.path.exists(GTDB_MERGED):
        with open(GTDB_MERGED, encoding="utf-8", errors="replace") as f:
            yield from f
        return
    for part in GTDB_PARTS:
        if not os.path.exists(part):
            continue
        with gzip.open(part, "rt", encoding="utf-8", errors="replace") as f:
            yield from f


def gtdb_table() -> dict[str, tuple[str, str]]:
    """``{numeric accession body: (lineage, species)}`` from GTDB r214."""
    global _gtdb_cache
    if _gtdb_cache is not None:
        return _gtdb_cache
    table: dict[str, tuple[str, str]] = {}
    for line in _iter_gtdb_lines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 2:
            continue
        key = normalise(parts[0])
        if not key or key in table:
            continue
        lineage = parts[1].strip()
        m = re.search(r"s__([^;]+)", lineage)
        table[key] = (lineage, m.group(1).strip() if m else "")
    _gtdb_cache = table
    return table


def _iter_ncbi_pairs():
    if os.path.exists(NCBI_SLIM):
        with open(NCBI_SLIM, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] and parts[1]:
                    yield parts[0], parts[1]
        return
    # Fall back to the raw NCBI report: col 1 = accession, col 8 = organism_name
    if os.path.exists(NCBI_BULK):
        with open(NCBI_BULK, encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 9 and parts[0] and parts[7]:
                    yield parts[0], parts[7]


def ncbi_table() -> dict[str, str]:
    """``{accession: organism_name}`` from the local NCBI RefSeq report."""
    global _ncbi_cache
    if _ncbi_cache is not None:
        return _ncbi_cache
    table: dict[str, str] = {}
    for acc, name in _iter_ncbi_pairs():
        table.setdefault(acc, name)
    # results previously fetched from the API, so a warm cache needs no network
    if os.path.exists(API_CACHE):
        with open(API_CACHE, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] and parts[1]:
                    table.setdefault(parts[0], parts[1])
    _ncbi_cache = table
    return table


def _ncbi_body_table() -> dict[str, str]:
    """``{numeric accession body: organism_name}`` -- version-insensitive index.

    Built lazily (it only costs anything when an accession misses the exact
    table) so a genome whose filename dropped the version still resolves.
    """
    global _ncbi_body_cache
    if _ncbi_body_cache is None:
        _ncbi_body_cache = {}
        for acc, name in ncbi_table().items():
            body = normalise(acc)
            if body:
                _ncbi_body_cache.setdefault(body, name)
    return _ncbi_body_cache


#: Second word of an NCBI organism name when it is a qualifier, not the species
#: epithet (``Escherichia coli str. K-12``, ``Bacillus sp. 1``).
_QUALIFIER = {"sp.", "sp", "spp.", "spp", "cf.", "cf", "str.", "str", "strain",
              "subsp.", "subsp", "ssp.", "bv.", "bv", "serovar", "pv."}


def species_name(name: str) -> str:
    """MIMEdb-friendly species-level name from an NCBI organism name.

    MIMEdb's ``species`` column holds bare binomials, while NCBI's
    ``organism_name`` appends the strain. Measured against
    ``database/mimedb/mimedb_microbes_v2.csv``: 'Clostridium perfringens ATCC
    13124' matches 0 rows, 'Clostridium perfringens' matches 4. So the strain
    qualifier has to come off before the lookup, not after.

        ``Clostridium perfringens ATCC 13124``        -> Clostridium perfringens
        ``Escherichia coli str. K-12 substr. MG1655`` -> Escherichia coli
        ``[Eubacterium] rectale ATCC 33656``          -> Eubacterium rectale
    """
    n = re.sub(r"\s+", " ", (name or "").replace("[", "").replace("]", "")).strip()
    words = n.split(" ")
    if not words or not words[0]:
        return ""
    if len(words) == 1:
        return words[0]
    if words[1].lower() in _QUALIFIER:
        return words[0]
    # ``Candidatus Xyz abc`` -- keep the placeholder plus the genus.
    if words[0].lower() == "candidatus" and len(words) >= 3:
        return f"{words[0]} {words[1]}"
    return f"{words[0]} {words[1]}"


def lookup_gtdb(acc: str) -> tuple[str, str]:
    key = normalise(acc)
    if not key:
        return "", ""
    return gtdb_table().get(key, ("", ""))


def lookup_ncbi(acc: str) -> str:
    """NCBI organism name for *acc*, tolerating a missing/extra version suffix."""
    table = ncbi_table()
    if acc in table:
        return table[acc]
    body = normalise(acc)
    if not body:
        return ""
    return _ncbi_body_table().get(body, "")


# ---------------------------------------------------------------------------
# Online fallbacks
# ---------------------------------------------------------------------------
def _api_fetch(acc: str) -> str:
    """NCBI Datasets API: accession -> organism_name. '' on any failure."""
    import json
    try:
        req = urllib.request.Request(
            DATASETS_API.format(acc=acc),
            headers={"User-Agent": "GEMAgent/2 taxonomy resolver"})
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = json.load(resp)
    except (urllib.error.URLError, ValueError, OSError, TimeoutError):
        return ""
    for report in data.get("reports") or []:
        name = (report.get("organism") or {}).get("organism_name")
        if name:
            return str(name).strip()
    return ""


def _api_cached(acc: str, allow_network: bool) -> str:
    """API lookup with an on-disk cache so reruns stay offline and fast."""
    from_cache = ""
    if os.path.exists(API_CACHE):
        with open(API_CACHE, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2 and parts[0] == acc:
                    from_cache = parts[1]
                    break
    if from_cache or not allow_network:
        return from_cache
    name = _api_fetch(acc)
    if name:
        os.makedirs(os.path.dirname(API_CACHE), exist_ok=True)
        with open(API_CACHE, "a", encoding="utf-8") as f:
            f.write(f"{acc}\t{name}\n")
        if _ncbi_cache is not None:
            _ncbi_cache.setdefault(acc, name)
    return name


def _gtdbtk_run(genome_file: str, output_dir: str) -> tuple[str, str]:
    """Last-resort GTDB-TK ``classify_wf`` for genomes absent from the tables.

    ``GTDBTK_DATA_PATH`` is only exported by the gtdbtk env's ``activate.d``,
    so the env has to be activated in a shell. Returns '' pairs on any failure
    -- a missing taxonomy must never abort a reconstruction.
    """
    import shutil
    import subprocess

    gtdbtk = shutil.which("gtdbtk") or paths._conda_env_bin("gtdbtk")
    if not gtdbtk or not genome_file:
        return "", ""
    work = os.path.join(output_dir, "gtdbtk")
    genome_dir = os.path.join(work, "input")
    summary = os.path.join(work, "output", "classify",
                           "gtdbtk.bac120.summary.tsv")
    alt = os.path.join(work, "output", "classify", "gtdbtk.ar53.summary.tsv")
    if not os.path.exists(summary) and not os.path.exists(alt):
        os.makedirs(genome_dir, exist_ok=True)
        dst = os.path.join(genome_dir, os.path.basename(genome_file))
        if not os.path.exists(dst):
            shutil.copy2(genome_file, dst)
        activate = os.path.join(os.path.dirname(os.path.dirname(gtdbtk)),
                                "..", "etc", "profile.d", "conda.sh")
        cmd = (f"source '{os.path.abspath(activate)}' && conda activate "
               f"'{os.path.basename(os.path.dirname(os.path.dirname(gtdbtk)))}' "
               f"&& '{gtdbtk}' classify_wf --genome_dir '{genome_dir}' "
               f"--out_dir '{os.path.join(work, 'output')}' "
               f"--extension {os.path.basename(genome_file).split('.', 1)[-1]} "
               f"--cpus 8")
        print(f"  [taxonomy] GTDB-TK run (this takes 10-30 min): {gtdbtk}")
        try:
            subprocess.run(["bash", "-lc", cmd], check=True,
                           stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
        except (subprocess.CalledProcessError, OSError) as e:
            print(f"  [taxonomy] GTDB-TK failed ({e}); continuing without.")
            return "", ""
    for path in (summary, alt):
        if not os.path.exists(path):
            continue
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = f.readlines()
        if len(lines) >= 2 and "\t" in lines[1]:
            lineage = lines[1].split("\t")[1].strip()
            m = re.search(r"s__([^;]+)", lineage)
            return lineage, (m.group(1).strip() if m else "")
    return "", ""


# ---------------------------------------------------------------------------
# Public entry points
# ---------------------------------------------------------------------------
def resolve(genome_id: str, genome_file: str = "", output_dir: str = "",
            org_name: str = "", taxonomy_tree: str = "",
            allow_network: bool = True, gtdb_fallback: bool = True,
            verbose: bool = True) -> TaxonomyResult:
    """Resolve both naming fields for one genome.

    Explicit *org_name* / *taxonomy_tree* arguments always win; the tables only
    fill what is missing. An explicit *org_name* also suppresses the NCBI
    lookup, but never the GTDB lineage lookup.
    """
    header = _fasta_first_header(genome_file)
    accession = extract_accession(genome_id, os.path.basename(genome_file or ""),
                                  header)
    res = TaxonomyResult(genome_id=genome_id, accession=accession)
    sources: list[str] = []

    if taxonomy_tree:
        res.gtdb_lineage = taxonomy_tree
        m = re.search(r"s__([^;]+)", taxonomy_tree)
        res.gtdb_species = m.group(1).strip() if m else ""
        sources.append("cli:taxonomy_tree")
    elif accession:
        res.gtdb_lineage, res.gtdb_species = lookup_gtdb(accession)
        if res.gtdb_lineage:
            sources.append("gtdb_r214")

    if org_name:
        res.org_name = org_name
        sources.append("cli:org_name")
    elif accession:
        res.org_name = lookup_ncbi(accession)
        if res.org_name:
            sources.append("ncbi_refseq")
        else:
            res.org_name = _api_cached(accession, allow_network)
            if res.org_name:
                sources.append("ncbi_api")

    # last resort for the lineage (never for org_name: a GTDB name is worse
    # than no name, see the module docstring)
    if not res.gtdb_lineage and gtdb_fallback:
        res.gtdb_lineage, res.gtdb_species = _gtdbtk_run(genome_file, output_dir)
        if res.gtdb_lineage:
            sources.append("gtdbtk")

    # last resort for org_name: GTDB species with its GTDB-only suffix removed
    # (``Clostridium_S acetobutylicum`` -> ``Clostridium acetobutylicum``).
    if not res.org_name and res.gtdb_species:
        stripped = re.sub(r"_(?:[A-Z]|[A-Z]{2,3}\d*)$", "", res.gtdb_species)
        stripped = re.sub(r"\s+sp\d+$", "", stripped).strip()
        if stripped:
            res.org_name = stripped
            sources.append("gtdb_species(stripped)")

    # MIMEdb keys on bare binomials, so a strain-qualified NCBI name would miss;
    # see species_name().
    res.species_name = species_name(res.org_name) or res.org_name

    res.source = "+".join(sources)
    if verbose:
        print("  [taxonomy] " + res.describe().replace("\n", "\n  [taxonomy] "))
        if not res.gtdb_lineage:
            print("  [taxonomy] WARNING: no GTDB lineage resolved "
                  "(reviews will run without taxonomic context)")
        if not res.org_name or res.org_name == genome_id:
            print("  [taxonomy] WARNING: no organism name resolved "
                  "(MIMEdb/literature context will be empty)")
    return res


def write_result(output_dir: str, res: TaxonomyResult) -> str:
    """Persist *res* to ``<output_dir>/taxonomy/taxonomy.tsv``."""
    out_dir = os.path.join(output_dir, OUT_SUBDIR)
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, OUT_FILE)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\t".join(OUT_COLUMNS) + "\n")
        f.write("\t".join(res.as_row()) + "\n")
    return path


def read_result(output_dir: str) -> TaxonomyResult | None:
    """Read back ``<output_dir>/taxonomy/taxonomy.tsv`` (None when absent)."""
    path = os.path.join(output_dir, OUT_SUBDIR, OUT_FILE)
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.read().splitlines()
    if len(lines) < 2:
        return None
    values = lines[1].split("\t")
    res = TaxonomyResult()
    for col, val in zip(OUT_COLUMNS, values):
        setattr(res, col, val)
    return res


def lineage_of(output_dir: str) -> str:
    """Convenience for the review steps: GTDB lineage for a finished resolve."""
    res = read_result(output_dir)
    return res.gtdb_lineage if res else ""


def org_name_of(output_dir: str) -> str:
    """Convenience for the review steps: NCBI organism name."""
    res = read_result(output_dir)
    return res.org_name if res else ""


def species_name_of(output_dir: str) -> str:
    """Convenience for MIMEdb lookups: species-level name."""
    res = read_result(output_dir)
    if not res:
        return ""
    return res.species_name or res.org_name


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description="Resolve taxonomy / organism name")
    ap.add_argument("--genome_id", default="")
    ap.add_argument("--genome_file", default="")
    ap.add_argument("--output_dir", default="")
    ap.add_argument("--org_name", default="")
    ap.add_argument("--taxonomy_tree", default="")
    ap.add_argument("--no_network", action="store_true")
    ap.add_argument("--no_gtdb_fallback", action="store_true")
    ap.add_argument("--write", action="store_true",
                    help="write <output_dir>/taxonomy/taxonomy.tsv")
    a = ap.parse_args()

    r = resolve(a.genome_id, a.genome_file, a.output_dir, a.org_name,
                a.taxonomy_tree, allow_network=not a.no_network,
                gtdb_fallback=not a.no_gtdb_fallback)
    if a.write and a.output_dir:
        print(f"\n  wrote {write_result(a.output_dir, r)}")
    print(f"\n  GTDB lineage : {r.gtdb_lineage or '<none>'}")
    print(f"  org_name     : {r.org_name or '<none>'}")
    print(f"  species_name : {r.species_name or '<none>'}")
