#!/usr/bin/env python3
"""
GEMAgent Ã¢â‚¬â€ Evidence + LLM/DR curated genome-scale model reconstruction.

What it does:
  genome FASTA Ã¢â€ â€™ gapseq evidence (BLAST/pathway/transporters) Ã¢â€ â€™ gapseq draft
  SBML Ã¢â€ â€™ 3 LLM reviews (biomass / gapfilling / GPR) Ã¢â€ â€™ curated GEM.

What it is for:
  Use gapseq's draft construction as the front-end, then inject traceable LLM
  curation at the four points where gapseq's automatic rules are known to be
  brittle. Gene essentiality quality is the primary proven benchmark.

Pipeline (native gapseq front-end + LLM reviews):
  0. Prompts (MIMEdb context) -- organism name + GTDB lineage resolved from the
                                    assembly accession by code/taxonomy.py
  1. BLAST + pathway topology        (gapseq find)
  2. Transporter prediction          (gapseq find-transport)
  3. Draft model                     (gapseq draft -- it also picks the biomass
                                      template, which is read back afterwards)
  4. KOFam annotation                (prodigal + KOFamScan; evidence for the
                                      biomass and GPR reviews)
  5. Growth medium prediction
  6. [Review 2] Biomass review       (LLM audits biomass precursors ->
                                      -draft-biomassadj.RDS)
  7. Gapfilling                      (code/gapfilling, a copy of gapseq's own
                                      gf.suite.R -- reads gapseq's `dat/` via the
                                      code/dat symlink and takes -c/-g weights)
  8. [Review 3] Gapfilling review    (LLM reviews the bypasses gapfill added, and
                                      re-runs gapfill with the updated excludes)
  9. [Review 3] GPR review           (LLM removes over-attached isozymes)
 10. [Review 1] Diffusion review     (NOT WIRED YET: gapseq injects its full
                                      diffusion sink set unconditionally, so this
                                      will prune the disapproved sinks after the
                                      filled model exists)

The active review stages are biomass, gapfilling, and GPR.

Resuming: --from_draft skips steps 1-3 when <id>-draft.RDS/.xml are already
present, so the reviews can be iterated on without redoing the ~4.5 h front end.
"""

import os
import sys
import re
import argparse
import shutil
import subprocess
import glob

# ==============================================================================
# Runtime bootstrap - shared library + external tool PATH resolution.
# ------------------------------------------------------------------------------
# (1) gemagent conda env ships a newer libstdc++ (CXXABI_1.3.15), but the process
#     otherwise loads the system /lib64/libstdc++.so.6 (too old), which breaks
#     cobra / sqlite3 / icu imports (used by the GPR-review scripts). Prepend the
#     conda env's lib dir to LD_LIBRARY_PATH so every child process inherits it.
# (2) blast/diamond live outside the conda env (Software_mdmclean); find.sh calls
#     them as bare commands, so prepend their dirs to PATH.
# ==============================================================================
_env_bin = os.path.dirname(os.path.abspath(sys.executable))
_env_root = os.path.dirname(_env_bin)
_env_lib = os.path.join(_env_root, "lib")
if os.path.isdir(_env_lib):
    _cur_ld = os.environ.get("LD_LIBRARY_PATH", "")
    if _env_lib not in _cur_ld.split(os.pathsep):
        os.environ["LD_LIBRARY_PATH"] = os.pathsep.join(
            p for p in (_env_lib, _cur_ld) if p
        )

_ext_bins = os.environ.get("GEMAGENT_EXTRA_BIN", "").split(os.pathsep)
for _b in _ext_bins:
    _b = _b.strip()
    _cur_path = os.environ.get("PATH", "")
    if _b and os.path.isdir(_b) and _b not in _cur_path.split(os.pathsep):
        os.environ["PATH"] = _b + os.pathsep + _cur_path

# ==============================================================================
# Global configuration
# ==============================================================================
n_threads = os.cpu_count() or 1
print("=" * 60)
print(f"  GEMAgent ({n_threads} threads)")
print("=" * 60)

VERBOSE       = 1
BITCUTOFF     = 200          # GapSeq-default discovery threshold
TAXONOMY      = "auto"
MIN_BS_CORE   = 100
FORCE_OFFLINE = False

# Review concurrency / retry knobs. Named here (env-overridable) instead of being
# buried as literals inside the subprocess argument lists, so the LLM load can be
# tuned per machine without editing code.
BIOMASS_WORKERS = int(os.environ.get("GEMAGENT_BIOMASS_WORKERS", "8"))
GPR_WORKERS     = int(os.environ.get("GEMAGENT_GPR_WORKERS", "16"))
REVIEW_RETRIES  = int(os.environ.get("GEMAGENT_REVIEW_RETRIES", "2"))

OUTPUT_SUBDIR = "gemagent_output"
_HERE         = os.path.dirname(os.path.abspath(__file__))
# gemagent.py may live at the project root (canonical) or inside code/.
CODE_DIR      = os.path.join(_HERE, "impl")

# Paths / external binaries all come from code/paths.py -- no absolute paths are
# hardcoded here. PYTHON defaults to the running interpreter, RSCRIPT to its
# sibling in the same conda env, and GAPSEQ_BIN is searched for (PATH, then any
# conda env) because gapseq lives in its own env.
sys.path.insert(0, CODE_DIR)
import paths  # noqa: E402
import taxonomy  # noqa: E402

PYTHON        = paths.PYTHON
RSCRIPT       = paths.RSCRIPT
# Native gapseq CLI. v2 deliberately does NOT vendor or patch gapseq's scripts:
# the front end (find / find-transport / draft) goes through this binary exactly
# as `gapseq doall` does, so the draft is the same class of model as the
# benchmark `models/gapseq/*.RDS`.
GAPSEQ_BIN    = paths.GAPSEQ_BIN
# `gapseq draft` thresholds Ã¢â‚¬â€ gapseq doall's own defaults.
DRAFT_HIGH_EVI_BS = os.environ.get("GEMAGENT_DRAFT_HIGH_EVI_BS", "200")
DRAFT_MIN_BS_CORE = os.environ.get("GEMAGENT_DRAFT_MIN_BS_CORE", "100")

# `gapseq find` thread count. find.sh splits each reaction query into 10 parts
# of ~0.6 KB and runs one BLAST per part; that is far below its `parallel --pipe
# --block 50k` chunk size, so its parallel branch degenerates to one job per part
# plus 0.3-2.6 s of spawn overhead. Measured with the env's BLAST 2.17.0:
#   serial tblastn             0.137-0.181 s/part   (1-seq)   2.445 s/part (90 seq)
#   parallel --pipe (default)  0.895-7.135 s/part   (1-seq)   2.753 s/part (90 seq)
# `-K 1` makes find.sh take the serial branch (identical output; BLAST's
# -num_threads does not change results). find-transport keeps its own -K because
# it runs a single blastp against uniprot_sprot, where threads do help.
DRAFT_FIND_THREADS = os.environ.get("GEMAGENT_FIND_THREADS", "1")

# --- gapseq subprocess environment -------------------------------------------
# The `gapseq` wrapper script does NOT export PATH, while gapseq_find.sh calls
# `tblastn`, `parallel`, `exonerate` (fastaindex / fastafetch) and `Rscript` by
# bare name. If the gapseq env's bin is not first on PATH, three things go
# silently wrong on this cluster (all three measured on csr32, 2026-09-17):
#
#   1. `tblastn` resolves to
#        Project/Software_mdmclean/ncbi-blast-2.13.0+/bin/tblastn
#      a build whose *startup* stalls ~3.12 s per invocation here -- even
#      `tblastn -version` takes real=3.123s with user=0.006s sys=0.010s, and an
#      empty query costs the same. With ~3 BLAST calls per reaction over ~12k
#      reactions that is ~31 h of pure stall (observed: 13.7 s/query, i.e. ~46 h
#      for bsub). The env's BLAST 2.17.0 starts in ~0.1 s.
#   2. `Rscript` is absent from PATH entirely (anaconda3/bin/Rscript does not
#      exist, rc=127), so `complex_detection.R` never runs. find.sh then does
#      `query=$(readlink -f subunit_tmp.fasta)` on a file that was never
#      created, and every query whose reference sequences mention
#      subunit/chain/polypeptide/component is silently mis-processed.
#      Symptom: `iterations>1` is 0 in subunits.log, while gapseq's own bsub
#      table has 678 of them.
#   3. `parallel` is absent, so find.sh takes its serial fallback branch
#      (`tblastn` without `-num_threads`) instead of `parallel --pipe -j N`.
#
# Prepending the gapseq env's bin fixes all three without touching gapseq.
GAPSEQ_ENV = dict(os.environ)
# Only prepend when GAPSEQ_BIN actually resolved: os.path.abspath("") is the cwd,
# and prepending that to PATH would silently shadow real tools.
_GAPSEQ_ENV_BIN = os.path.dirname(os.path.abspath(GAPSEQ_BIN)) if GAPSEQ_BIN else ""
_GAPSEQ_ENV_PREFIX = _GAPSEQ_ENV_BIN + os.sep if _GAPSEQ_ENV_BIN else ""
if _GAPSEQ_ENV_BIN and not GAPSEQ_ENV.get("PATH", "").startswith(_GAPSEQ_ENV_PREFIX):
    GAPSEQ_ENV["PATH"] = _GAPSEQ_ENV_BIN + os.pathsep + GAPSEQ_ENV.get("PATH", "")


def report_gapseq_tools():
    """Print which tblastn / Rscript / parallel / exonerate gapseq will use.

    Cheap insurance against the silent-degradation mode described above: if
    `tblastn` is not the gapseq env's BLAST, find.sh will be ~30x slower on
    startup alone.
    """
    if not GAPSEQ_BIN:
        print("\n  [Warning] gapseq CLI not found; the gapseq front end "
              "(find / find-transport / draft) will not run.\n"
              "            Set GAPSEQ_BIN or put `gapseq` on PATH.")
        return
    print(f"\n  gapseq CLI: {GAPSEQ_BIN}")
    print("  gapseq tool resolution (env bin first on PATH):")
    for tool in ("tblastn", "blastp", "makeblastdb", "Rscript",
                 "parallel", "fastaindex", "fastafetch", "exonerate"):
        found = shutil.which(tool, path=GAPSEQ_ENV["PATH"])
        mark = "ok " if found and found.startswith(_GAPSEQ_ENV_PREFIX) else "!! "
        print(f"    {mark}{tool:<12} {found or '<NOT FOUND>'}")

PROJECT_ROOT  = paths.PROJECT_ROOT

# --- gapfill medium policy ----------------------------------------------------
# Which medium should gapfilling run on?
#
# The benchmark models (gapseqEval's models/gapseq/*.RDS) were NOT gapfilled on
# gapseq's complete medium (database/media/ALLmed.csv). That was tested and
# disproved: rerunning `gapseq doall -m ALLmed.csv` for mgen still yields a model
# that cannot grow on Complete2 (growth = 0). The one medium that reproduces the
# benchmark model's behaviour is the organism's *evaluation* medium:
#
#   gapfill medium                 mgen growth on Complete2
#   ----------------------------   -----------------------
#   predicted (doall default)      0.000000
#   ALLmed.csv ("complete")        0.000000
#   Complete2 (the eval medium)    0.398295   <- reproduces the benchmark
#
# Cross-check: all five benchmark RDS models cover 100% of their own evaluation
# medium (ecol GS_MM_glc, paer GS_MM_succ, bsub/sone LB_marinos, mgen Complete2).
#
# Why it matters: gapfilling only guarantees growth on the medium it is run on.
# If that medium differs from the evaluation medium, a model can be degenerate
# there (growth = 0) -> every knockout looks lethal -> all genes "essential" ->
# a meaningless F1. mgen is the sharp case: its predicted medium supplies
# pantothenate + cysteine (so gapfill builds "pantothenate + Cys -> CoA"), while
# Complete2 supplies CoA + Cys-Gly instead (so the built route is dead). Running
# gapfill on the evaluation medium makes gapfill add the CoA/Cys-Gly uptake
# directly (as the benchmark model does).
#
#   eval      -> use the medium passed via --eval_medium (default when given)
#   predicted -> <output_dir>/<genome_id>-medium.csv   (legacy GEMAgent behaviour)
#   complete  -> database/media/ALLmed.csv             (kept for experiments)
GAPFILL_MEDIUM_MODE = os.environ.get(
    "GEMAGENT_GAPFILL_MEDIUM", ""
).strip().lower()
GAPFILL_COMPLETE_MEDIUM = paths.db("media", "ALLmed.csv")


# ==============================================================================
# Helpers
# ==============================================================================

def run_command(cmd, description="Running command", check=True, env=None):
    print(f"\n{'Ã¢â€â‚¬' * 50}")
    print(f"  {description}")
    print(f"  CMD: {' '.join(cmd) if isinstance(cmd, list) else cmd}")
    print(f"{'Ã¢â€â‚¬' * 50}")
    try:
        if isinstance(cmd, str):
            result = subprocess.run(cmd, shell=True, check=check,
                                    capture_output=True, text=True, env=env)
        else:
            result = subprocess.run(cmd, check=check,
                                    capture_output=True, text=True, env=env)
        if result.stdout:
            print(result.stdout)
        if result.stderr:
            print("  STDERR:", result.stderr[:500], file=sys.stderr)
        return result
    except subprocess.CalledProcessError as e:
        print(f"  Error: {e}", file=sys.stderr)
        if e.stdout:
            print(f"  STDOUT: {e.stdout[:500]}", file=sys.stderr)
        if e.stderr:
            print(f"  STDERR: {e.stderr[:500]}", file=sys.stderr)
        if check:
            raise
        return None


def find_genome_file(input_dir):
    for ext in ['*.fna', '*.fasta', '*.fa', '*.fna.gz', '*.fasta.gz', '*.fa.gz']:
        files = glob.glob(os.path.join(input_dir, ext))
        if files:
            return files[0]
    for ext in ['*.faa', '*.faa.gz']:
        files = glob.glob(os.path.join(input_dir, ext))
        if files:
            print(f"  Warning: found protein file {files[0]}")
            return files[0]
    raise FileNotFoundError(f"No genome file found in {input_dir}")


def extract_genome_id(filepath):
    name = os.path.basename(filepath)
    for ext in ['.fna.gz', '.fasta.gz', '.fa.gz', '.fna', '.fasta', '.fa',
                '.faa.gz', '.faa']:
        if name.endswith(ext):
            name = name[:-len(ext)]
            break
    return name


def _is_stale(output_file, input_files):
    """Return True when output is missing or older than any existing input."""
    if not os.path.exists(output_file):
        return True
    output_mtime = os.path.getmtime(output_file)
    for path in input_files:
        if path and os.path.exists(path) and os.path.getmtime(path) > output_mtime:
            return True
    return False


# ==============================================================================
# Cross-tool seeding
# ==============================================================================

def seed_draft_from_sbml(output_dir, genome_id, sbml_path):
    """Seed <output_dir>/<genome_id>-draft.{RDS,xml} from an external SBML.

    Cross-tool entry point (`--draft_sbml`): lets GEMAgent graft its LLM reviews
    onto a model built by another tool (ModelSEED, CarveMe, ...) instead of a
    gapseq draft. The SBML is expected to carry gapseq-style coordinate gene
    names already (see cross_tool/relabel_model.py); the RDS is materialised by
    seed_foreign_draft.R, which also undoes the SBML-id escape that would
    otherwise destroy the contig/coordinate structure that the GPR review's
    KOFam bridge and the official gapseqEval parse out of gene names.

    Returns True when both the RDS and the XML are in place.
    """
    draft_rds = os.path.join(output_dir, f"{genome_id}-draft.RDS")
    draft_xml = os.path.join(output_dir, f"{genome_id}-draft.xml")
    if not os.path.isfile(sbml_path):
        print(f"  Error: --draft_sbml file not found: {sbml_path}", file=sys.stderr)
        return False
    if (os.path.exists(draft_rds) and os.path.exists(draft_xml) and
            os.path.getmtime(draft_rds) >= os.path.getmtime(sbml_path)):
        print("  Draft already seeded from the external SBML Ã¢â‚¬â€ skipping.")
        return True
    print(f"  Seeding draft from external SBML: {sbml_path}")
    res = run_command([
        RSCRIPT, f"{CODE_DIR}/agent/agent_gpr_check/seed_foreign_draft.R",
        "-x", sbml_path, "-o", draft_rds, "-v",
    ], "Step 0: Seed draft from external SBML", check=False)
    if (res is None or res.returncode != 0 or
            not os.path.isfile(draft_rds) or os.path.getsize(draft_rds) == 0):
        print("  Seeding the draft RDS from the external SBML failed.", file=sys.stderr)
        return False
    # The XML is the same model; keep a copy under the name the pipeline expects
    # (biomass evidence, GPR prompts and the SBML-only stages read <id>-draft.xml).
    shutil.copy2(sbml_path, draft_xml)
    print(f"  Draft seeded: {os.path.basename(draft_rds)} "
          f"+ {os.path.basename(draft_xml)}")
    return True


# ==============================================================================
# Pipeline steps
# ==============================================================================

def step_prompts(output_dir, org_name):
    run_command([
        PYTHON, f"{CODE_DIR}/prompt/org_prompt_gen.py",
        "--output_dir", output_dir, "--org_name", org_name,
    ], "Step 0: Generate prompts")


def step_find_pathways(genome_file, output_dir, genome_id):
    pathway_file = os.path.join(output_dir, f"{genome_id}-all-Pathways.tbl")
    reaction_file = os.path.join(output_dir, f"{genome_id}-all-Reactions.tbl")
    if os.path.exists(pathway_file) and os.path.exists(reaction_file):
        print("  Pathway files already exist Ã¢â‚¬â€ skipping.")
        return True
    cmd = [
        GAPSEQ_BIN, "find",
        "-f", output_dir + "/",
        "-b", str(BITCUTOFF), "-p", "all", "-t", TAXONOMY,
        # serial branch on purpose -- see DRAFT_FIND_THREADS
        "-K", str(DRAFT_FIND_THREADS), "-v", str(VERBOSE),
    ]
    if FORCE_OFFLINE:
        cmd.append("-O")
    cmd.append(genome_file)
    report_gapseq_tools()
    run_command(cmd, "Step 1: Find pathways (gapseq find)",
                env=GAPSEQ_ENV)
    return os.path.exists(pathway_file) and os.path.exists(reaction_file)


def step_predict_transporters(genome_file, output_dir, genome_id):
    transporter_file = os.path.join(output_dir, f"{genome_id}-Transporter.tbl")
    if os.path.exists(transporter_file):
        print("  Transporter file already exists Ã¢â‚¬â€ skipping.")
        return True
    run_command([
        GAPSEQ_BIN, "find-transport",
        "-f", output_dir + "/",
        "-b", str(BITCUTOFF), "-K", str(n_threads), "-v", str(VERBOSE),
        genome_file,
    ], "Step 2: Predict transporters (gapseq find-transport)",
       env=GAPSEQ_ENV)
    return os.path.exists(transporter_file)


# The biomass template gapseq picked is stamped verbatim into the draft's biomass
# reaction name, so the finished draft -- not a separate 16S or MIMEdb call -- is
# the authoritative source for it. gapseq's own `-b auto` logic (16S classifier,
# falling back to a k-mer classifier) already decided it.
_DRAFT_BIOMASS_TEMPLATE = {
    "bacterialgrampositivebiomassreaction": "Gram_pos",
    "bacterialgramnegativebiomassreaction": "Gram_neg",
    "archaealbiomassreaction": "archaea",
}


def step_read_biomass_template(output_dir, genome_id):
    """Recover the biomass template from the draft gapseq built.

    gapseq writes the choice into the biomass reaction name, e.g. "Bacterial
    Gram-negative biomass reaction" for biomass_Gram_neg.json. Reading it back
    keeps GEMAgent's downstream biomass review aligned with what gapseq actually
    used, instead of re-deriving the call from a separate classifier.
    """
    xml = os.path.join(output_dir, f"{genome_id}-draft.xml")
    if not os.path.exists(xml):
        print(f"  Warning: draft SBML not found: {xml}")
        return ""
    # gapseq writes the biomass reaction *after* the whole <species> block, i.e.
    # ~7.7 MB into a 10.4 MB document, so this must not be a fixed-size head
    # read. A single line scan is correct and cheap.
    name = ""
    with open(xml, encoding="utf-8", errors="replace") as f:
        for line in f:
            if "<reaction" not in line or "biomass reaction" not in line:
                continue
            m = re.search(r'\bname="([^"]*[Bb]iomass reaction)"', line)
            if m:
                name = m.group(1)
                break
    if not name:
        print("  Warning: no biomass reaction found in the draft")
        return ""
    tmpl = _DRAFT_BIOMASS_TEMPLATE.get(re.sub(r"[\s\-_]+", "", name.lower()), "")
    print(f"  Draft biomass reaction: {name!r} -> template {tmpl or 'UNKNOWN'}")
    return tmpl


# ---------------------------------------------------------------------------
# Ã¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢Â
# v2 core
# Ã¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢Â

def _biomass_stem(gram_staining):
    biomass_map = {
        "positive": "Gram_pos", "gram_pos": "Gram_pos", "pos": "Gram_pos",
        "gram_positive": "Gram_pos",
        "negative": "Gram_neg", "gram_neg": "Gram_neg", "neg": "Gram_neg",
        "gram_negative": "Gram_neg",
        "archaea": "archaea",
    }
    return biomass_map.get(str(gram_staining).strip().lower(), "Gram_pos")


def _biomass_template(gram_staining):
    """gapseq biomass template JSON for the given gram-staining result."""
    return paths.db("biomass", f"biomass_{_biomass_stem(gram_staining)}.json")


def step_build_draft_gapseq(output_dir, genome_id, genome_file, force=False):
    """Build the draft model with `gapseq draft` (no vendored script, no patch).

    gapseq's own draft step selects reactions from BLAST + pathway-topology
    evidence, exactly as `gapseq doall` does, so the result is the same class of
    model as the benchmark `models/gapseq/*.RDS`. `-b auto` lets gapseq choose
    the biomass template; step_read_biomass_template() reads that choice back out
    of the draft afterwards.
    """
    rds_file = os.path.join(output_dir, f"{genome_id}-draft.RDS")
    if os.path.exists(rds_file) and not force:
        print("  Draft model already exists Ã¢â‚¬â€ skipping.")
        return True

    blast_file = os.path.join(output_dir, f"{genome_id}-all-Reactions.tbl")
    trans_file = os.path.join(output_dir, f"{genome_id}-Transporter.tbl")
    pwy_file   = os.path.join(output_dir, f"{genome_id}-all-Pathways.tbl")
    if not os.path.exists(blast_file):
        print(f"  Error: BLAST file not found: {blast_file}")
        return False

    run_command([
        GAPSEQ_BIN, "draft",
        "-r", blast_file,
        "-t", trans_file if os.path.exists(trans_file) else os.devnull,
        "-p", pwy_file if os.path.exists(pwy_file) else os.devnull,
        "-c", genome_file,
        "-u", str(DRAFT_HIGH_EVI_BS), "-l", str(DRAFT_MIN_BS_CORE),
        "-b", "auto",
        "-n", genome_id,
        "-f", output_dir + "/",
    ], "Step 4: Build draft model (gapseq draft)",
       env=GAPSEQ_ENV)
    return os.path.exists(rds_file)


#: Files the GPR review needs alongside the KO table to turn it into per-gene
#: evidence. The tsv alone is useless on its own: it is keyed by prodigal
#: ordinal, and mapping that back to the model's coordinate-style gene ids needs
#: either the prodigal GFF or the protein FASTA (whose headers carry the
#: coordinates). A directory holding only the tsv -- e.g. one reused from an
#: older run -- used to be accepted as "already annotated", which silently
#: produced a 0%-coverage prompt and an all-UNCERTAIN GPR review.
_KO_COORD_SOURCES = ("_proteins.faa", "_prodigal.gff")


def _kofam_coord_source(kofam_tsv):
    """Sibling faa/gff that provides gene coordinates, or '' when absent."""
    base = os.path.basename(kofam_tsv)
    if not base.endswith("_kofam.tsv"):
        return ""
    stem = base[: -len("_kofam.tsv")]
    for suffix in _KO_COORD_SOURCES:
        cand = os.path.join(os.path.dirname(kofam_tsv), stem + suffix)
        if os.path.exists(cand) and os.path.getsize(cand) > 0:
            return cand
    return ""


def _find_kofam_tsv(output_dir, genome_id):
    """Ã¦Å¸Â¥Ã¦â€°Â¾ KOFamScan Ã¦Â³Â¨Ã©â€¡Å  tsvÃ¯Â¼Å’Ã¥â€¦Â¼Ã¥Â®Â¹Ã¥Â¤Å¡Ã§Â§ÂÃ§â€ºÂ®Ã¥Â½â€¢/Ã¥â€˜Â½Ã¥ÂÂÃ§ÂºÂ¦Ã¥Â®Å¡Ã¯Â¼Ë†Ã¦Å“ÂªÃ¦â€°Â¾Ã¥Ë†Â°Ã¨Â¿â€Ã¥â€ºÅ¾Ã§Â©ÂºÃ¤Â¸Â²Ã¯Â¼â€°Ã£â‚¬â€š"""
    input_dir = os.path.dirname(output_dir.rstrip("/"))
    patterns = [
        os.path.join(output_dir, "kofam", f"{genome_id}_kofam.tsv"),
        os.path.join(output_dir, f"{genome_id}_kofam.tsv"),
        os.path.join(input_dir, "*_kofam", f"{genome_id}_kofam.tsv"),
        os.path.join(input_dir, f"{genome_id}_kofam", f"{genome_id}_kofam.tsv"),
    ]
    for pat in patterns:
        hits = glob.glob(pat)
        if hits:
            return hits[0]
    return ""


def step_annotate_kofam(output_dir, genome_id, genome_file, force=False):
    """Gene functional annotation (prodigal + KOFamScan).

    Provides gene-level KO evidence for the GPR review (step_gpr_check).
    Idempotent: skips when the annotation is complete, i.e. the KO table *and*
    its gene-coordinate companion (faa/gff) are both present and non-empty.
    """
    existing = _find_kofam_tsv(output_dir, genome_id)
    if existing and not force:
        if os.path.getsize(existing) > 0 and _kofam_coord_source(existing):
            print(f"  KOFam annotation already present: {existing}")
            return existing
        why = ("it is empty" if os.path.getsize(existing) == 0
               else "its gene-coordinate faa/gff companion is missing")
        print(f"  KOFam table found but unusable ({existing}): {why}. "
              f"Re-running annotation so the GPR review has KO evidence.")

    kofam_dir = os.path.join(output_dir, "kofam")
    os.makedirs(kofam_dir, exist_ok=True)
    target = os.path.join(kofam_dir, f"{genome_id}_kofam.tsv")
    script = os.path.join(CODE_DIR, "agent", "agent_gpr_check", "annotate_kofam.sh")

    if not os.path.isfile(script):
        print(f"  [Warning] annotate_kofam.sh not found: {script}")
        return ""

    # Fail loudly on a missing toolchain *before* spending time on prodigal: the
    # GPR review is blind without KO evidence (that is exactly how its KO
    # coverage once silently dropped to 0%).
    absent = paths.missing_kofam_tools()
    if absent:
        print("  [Warning] KOFam annotation skipped; missing:")
        for item in absent:
            print(f"            {item}")
        print("            Set GEMAGENT_KOFAM_DIR (or KOFAM_PROFILE / KOFAM_KO_LIST) "
              "and/or put prodigal + exec_annotation on PATH.")
        return ""

    tool_env = dict(os.environ)
    tool_env.update({k: v for k, v in paths.KOFAM_ENV.items() if v})
    run_command(
        ["bash", script, genome_file, kofam_dir, genome_id],
        "Step 5k: Annotate genes (prodigal + KOFamScan)",
        check=False,
        env=tool_env,
    )
    if os.path.exists(target) and os.path.getsize(target) > 0:
        print(f"  KOFam annotation done: {target}")
        return target
    print("  [Warning] KOFam annotation failed; continuing without KO evidence.")
    return ""


def _gapfill_model_file(output_dir, genome_id):
    """Gapfill input = the biomass-adjusted draft.

    step_biomass_check rewrites <id>-draft.RDS *in place* and keeps the raw
    gapseq draft as <id>-draft-nobiomass.RDS.  This matters because gf.suite.R
    derives its output id with gsub('-draft','',basename(model)) Ã¢â‚¬â€ so a model
    named <id>-draft.RDS yields the expected <id>-gapfilled.xml, whereas a
    differently-named copy (e.g. <id>-draft-biomassadj.RDS) would yield
    <id>-biomassadj-gapfilled.xml and break the downstream path contract.
    """
    raw   = os.path.join(output_dir, f"{genome_id}-draft-nobiomass.RDS")
    adj   = os.path.join(output_dir, f"{genome_id}-draft-biomassadj.RDS")
    draft = os.path.join(output_dir, f"{genome_id}-draft.RDS")
    if os.path.exists(raw):
        return draft            # adjustment was applied in place
    if os.path.exists(adj):
        return adj              # legacy layout
    return draft


def step_biomass_check(output_dir, genome_id, gram_staining, org_name="",
                       taxonomy_tree="", force=False, medium_csv=None,
                       kofam_tsv=None):
    """[Review 2] LLM biomass audit.

    1. biomass_evidence.py  Ã¢â‚¬â€ producibility gap per precursor (draft model + medium)
    2. biomass_check_prompt.py Ã¢â‚¬â€ pathway-centric prompts per template group
    3. review_biomass_check.py Ã¢â‚¬â€ LLM verdict KEEP / DROP / UNCERTAIN
    4. apply_biomass_decisions.py Ã¢â‚¬â€ writes <id>-draft-biomassadj.RDS (mass-renormalized)

    Returns True when a biomass-adjusted RDS was produced.
    """
    draft_xml = os.path.join(output_dir, f"{genome_id}-draft.xml")
    draft_rds = os.path.join(output_dir, f"{genome_id}-draft.RDS")
    raw_rds   = os.path.join(output_dir, f"{genome_id}-draft-nobiomass.RDS")
    raw_xml   = os.path.join(output_dir, f"{genome_id}-draft-nobiomass.xml")
    adj_rds   = os.path.join(output_dir, f"{genome_id}-draft-biomassadj.RDS")
    # The medium the audit reasons about must be the one gapfilling will run
    # on. <genome_id>-medium.csv is only produced by step_predict_medium(),
    # which is skipped when the caller passes --gapfill_medium/--eval_medium
    # explicitly -- so relying on it alone made biomass_evidence.py run with an
    # empty --medium_csv (0 compounds, wt_growth 0) and the LLM audit blind
    # (observed on cper 2026-09-20). Prefer the resolved medium, then fall back.
    medium    = medium_csv or os.path.join(output_dir, f"{genome_id}-medium.csv")
    if not medium or not os.path.exists(medium):
        print(f"  Warning: biomass evidence has no medium ({medium!r} missing); "
              f"producibility gaps will be judged without one.")
    bdir      = os.path.join(output_dir, "biomass_check")
    biomass_json = _biomass_template(gram_staining)

    if not os.path.exists(draft_xml) or not os.path.exists(draft_rds):
        print("  Skipping biomass review (no draft xml/RDS).")
        return False
    if os.path.exists(raw_rds) and not force:
        print("  Biomass-adjusted draft already applied in place Ã¢â‚¬â€ skipping.")
        return True
    os.makedirs(bdir, exist_ok=True)

    evidence = os.path.join(bdir, "biomass_evidence.json")
    run_command([
        PYTHON, f"{CODE_DIR}/agent/agent_biomass_check/biomass_evidence.py",
        "--model_xml", draft_xml,
        "--biomass_json", biomass_json,
        "--medium_csv", medium if os.path.exists(medium) else "",
        "--out_dir", bdir,
        "--kofam", kofam_tsv or _find_kofam_tsv(output_dir, genome_id),
        "--template_dir", os.path.dirname(biomass_json),
    ], "Step 7a: Biomass evidence (producibility gaps)", check=False)
    if not os.path.exists(evidence):
        print("  Biomass evidence failed Ã¢â‚¬â€ keeping original biomass.")
        return False

    run_command([
        PYTHON, f"{CODE_DIR}/agent/agent_biomass_check/biomass_check_prompt.py",
        "--evidence", evidence,
        "--output_dir", output_dir,
        "--out_dir", bdir,
        "--org_name", org_name or genome_id,
        "--taxonomy_tree", taxonomy_tree or "",
        "--audit_only",
    ], "Step 7b: Build biomass review prompts", check=False)

    run_command([
        PYTHON, f"{CODE_DIR}/agent/agent_biomass_check/review_biomass_check.py",
        "--prompt_dir", os.path.join(bdir, "prompt"),
        "--out_dir", bdir,
        "--evidence", evidence,
        "--workers", str(BIOMASS_WORKERS),
    ], "Step 7c: LLM biomass review", check=False)

    decisions = os.path.join(bdir, "biomass_decisions.json")
    if not os.path.exists(decisions):
        print("  Biomass decisions missing Ã¢â‚¬â€ keeping original biomass.")
        return False

    run_command([
        PYTHON, f"{CODE_DIR}/agent/agent_biomass_check/apply_biomass_decisions.py",
        "--draft_rds", draft_rds,
        "--biomass_json", biomass_json,
        "--decisions", decisions,
        "--out_dir", bdir,
        "--out_rds", adj_rds,
    ], "Step 7d: Apply biomass decisions", check=False)
    if not os.path.exists(adj_rds):
        print("  Biomass apply failed Ã¢â‚¬â€ keeping original biomass.")
        return False

    # Keep the raw gapseq draft, then rewrite the canonical draft *in place* so
    # that gf.suite.R emits <id>-gapfilled.* (see _gapfill_model_file).
    if not os.path.exists(raw_rds):
        shutil.copy2(draft_rds, raw_rds)
        if os.path.exists(draft_xml):
            shutil.copy2(draft_xml, raw_xml)
    shutil.move(adj_rds, draft_rds)
    for cand in (os.path.join(output_dir, f"{genome_id}-draft-biomassadj.xml"),
                 os.path.join(output_dir, f"{genome_id}-draft-biomassadj.RDS.xml")):
        if os.path.exists(cand):
            shutil.move(cand, draft_xml)
            break
    print(f"  Biomass-adjusted draft written in place: {os.path.basename(draft_rds)} "
          f"(raw kept as {os.path.basename(raw_rds)})")
    return True


# Ã¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢Â
# Post-construction
# Ã¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢ÂÃ¢â€¢Â

def step_predict_medium(output_dir, genome_id):
    medium_file = os.path.join(output_dir, f"{genome_id}-medium.csv")
    model_file = os.path.join(output_dir, f"{genome_id}-draft.RDS")
    pathway_file = os.path.join(output_dir, f"{genome_id}-all-Pathways.tbl")
    if not os.path.exists(model_file):
        print("  Warning: no draft model for medium prediction")
        return False
    if not _is_stale(medium_file, [model_file, pathway_file]):
        print("  Medium file is current Ã¢â‚¬â€ skipping.")
        return True
    run_command([
        GAPSEQ_BIN, "medium", "-m", model_file,
        "-p", pathway_file if os.path.exists(pathway_file) else os.devnull,
        "-f", output_dir + "/", "-o", os.path.basename(medium_file),
    ], "Step 10: Predict growth medium (official gapseq)", env=GAPSEQ_ENV)
    return os.path.exists(medium_file)


def _gapfill_exclude_file(output_dir, genome_id):
    """Path of the gapfill exclusion list, created empty if absent.

    ``gf.suite.R -d`` takes a TSV of reaction seeds to keep OUT of the gapfill
    candidate pool. The list is owned by the post-gapfill LLM review
    (``review_gapfill_added.py``), which appends its high-confidence bypass
    EXCLUDEs to this very file and then re-runs gapfilling; it reverts its
    appends if that re-run fails.

    The file is deliberately rebuilt empty at the start of every pipeline run
    rather than inheriting the previous run's contents: those came from a
    different LLM session (possibly a different thinking mode) and reusing them
    shrinks the candidate pool in ways the current run never justified, which
    can make gapfilling infeasible. The review regenerates its own exclusions
    within the run, so nothing is lost.

    (v1 seeded this list from ``annotated_candidates.json``, the output of a
    candidate-review stage that v2 no longer runs -- gemagent-v2's front end is
    plain ``gapseq find/find-transport/draft``. Reading that file was therefore
    removed along with the stage.)
    """
    exclude_file = os.path.join(output_dir, f"{genome_id}-gapfill-exclude.tsv")
    content = "seed\n"

    # Only rewrite when content actually changes, so the mtime-based staleness
    # check in step_gapfilling doesn't force gapfilling to re-run every time.
    old = ""
    if os.path.exists(exclude_file):
        with open(exclude_file) as f:
            old = f.read()
    if old != content:
        with open(exclude_file, "w") as f:
            f.write(content)
    return exclude_file


# NOTE: the v1 "gapfill failure diagnosis" step (iterative ddmin over the exclude
# list, code/diagnostics/find_minimal_blocking_excludes.R) was removed outright:
# that script never existed in v1 or v2, so the switch only ever printed a
# warning. Re-add the script before restoring the step.


def _resolve_gapfill_medium(output_dir, genome_id, gapfill_medium=None,
                            eval_medium=None):
    """Choose the medium gapfilling (and its LLM review) should run on.

    Precedence:
      1. an explicit ``gapfill_medium`` path always wins;
      2. otherwise an explicit ``GAPFILL_MEDIUM_MODE``;
      3. otherwise -- if an ``eval_medium`` (the medium this model will be
         evaluated on) is known, gapfill on it. This is what the benchmark
         models did, and it is the only choice that guarantees the model is not
         degenerate on the medium its F1 is measured on;
      4. otherwise fall back to the organism's predicted medium.
    """
    if gapfill_medium:
        return gapfill_medium
    mode = GAPFILL_MEDIUM_MODE
    if mode == "complete":
        return GAPFILL_COMPLETE_MEDIUM
    if mode == "predicted":
        return os.path.join(output_dir, f"{genome_id}-medium.csv")
    if mode == "eval":
        if not eval_medium:
            raise SystemExit("  Error: GAPFILL_MEDIUM_MODE=eval requires "
                             "--eval_medium <csv>")
        return eval_medium
    # unset -> prefer the evaluation medium when we know it
    if eval_medium:
        return eval_medium
    return os.path.join(output_dir, f"{genome_id}-medium.csv")


def step_gapfilling(output_dir, genome_id, gapfill_medium=None, model_rds=None,
                    eval_medium=None):
    gapfilled_xml = os.path.join(output_dir, f"{genome_id}-gapfilled.xml")
    model_file = model_rds or _gapfill_model_file(output_dir, genome_id)
    draft_xml = os.path.join(output_dir, f"{genome_id}-draft.xml")
    medium_file = _resolve_gapfill_medium(output_dir, genome_id, gapfill_medium,
                                          eval_medium)
    if not os.path.exists(model_file):
        print("  Error: draft model not found for gapfilling")
        return False

    exclude_file = _gapfill_exclude_file(output_dir, genome_id)

    # gapseq's gf.suite.R takes the reaction weights (-c) and the reaction->gene
    # table (-g) as explicit files. Both are produced by native `gapseq draft`.
    weights_file = os.path.join(output_dir, f"{genome_id}-rxnWeights.RDS")
    rxg_file = os.path.join(output_dir, f"{genome_id}-rxnXgenes.RDS")
    if not os.path.exists(weights_file) or not os.path.exists(rxg_file):
        print(f"  Error: gapfill weights missing ({weights_file} / {rxg_file})")
        return False

    if not _is_stale(gapfilled_xml, [model_file, draft_xml, medium_file,
                                     exclude_file, weights_file, rxg_file]):
        print("  Gapfilled model is current Ã¢â‚¬â€ skipping.")
        return True

    def _run_gapfill(label):
        # Official gapseq fill is invoked by the adapter. Exclusions are
        # applied by filtering the candidate tables, not by patching gapseq.
        from gemagent.adapters.gapseq import fill_with_excludes
        print(f"  {label}")
        return fill_with_excludes(
            output_dir, genome_id, model_file, medium_file, exclude_file,
            bcore=MIN_BS_CORE,
        )

    def _gapfill_succeeded(res):
        return bool(res) and os.path.exists(gapfilled_xml)

    result = _run_gapfill("Step 11: Gapfilling")
    if _gapfill_succeeded(result):
        return True

    # ---- gapfilling failed -------------------------------------------------
    # The first gapfill always runs with an empty exclusion list (see
    # _gapfill_exclude_file), so LLM exclusions cannot be the cause here: the
    # only inputs are the draft, the medium, the weights and the thresholds.
    # The most common cause is a medium that is too poor for the draft, so point
    # at that (--eval_medium aligns gapfilling with the medium the model will be
    # scored on).
    print("  Error: gapfilling failed (no feasible solution).")
    print("  The exclusion list was empty, so no LLM pruning is involved. Check "
          "the medium gapfilling ran on (pass --eval_medium to align it with the "
          "evaluation medium) and see the gf.suite.R output above.")
    return False


def step_review_gapfill_added(output_dir, genome_id, org_name="",
                              taxonomy_tree="", gapfill_medium=None,
                              max_iter=3, eval_medium=None,
                              retry_one_by_one=False):
    """Post-gapfill LLM bypass review (iterative).

    Reviews reactions gapfilling added (that were never candidate-reviewed),
    EXCLUDEs wrong-pathway/bypass leaks, and re-runs gapfilling with the updated
    exclusion list until the LLM stops excluding reactions.
    """
    model_file = _gapfill_model_file(output_dir, genome_id)
    medium_file = _resolve_gapfill_medium(output_dir, genome_id, gapfill_medium,
                                          eval_medium)
    rga_cmd = [
        PYTHON, f"{CODE_DIR}/agent/review_gapfill_added.py",
        "--output_dir", output_dir,
        "--genome_id", genome_id,
        "--org_name", org_name or genome_id,
        "--taxonomy_tree", taxonomy_tree or "",
        "--draft_model", model_file,
        "--medium", medium_file,
        "--bcore", str(MIN_BS_CORE),
        "--max_iter", str(max_iter),
    ]
    if retry_one_by_one:
        rga_cmd.append("--retry-one-by-one")
    result = run_command(rga_cmd, "Step 11b: Post-gapfill LLM bypass review",
                         check=False)
    final_decisions = os.path.join(output_dir, f"{genome_id}-gapfill-decisions-final.json")
    return result is not None and result.returncode == 0 and os.path.exists(final_decisions)


def _gpr_applier(model_path):
    """Pick the GPR applier from how the model stores its associations.

    fbc-style SBML (gapseq, AGORA) carries ``geneProductAssociation`` elements
    and is edited through cobra; CarveMe-style notes SBML keeps
    ``GENE_ASSOCIATION:`` inside ``<notes>`` and must be edited as text, because
    a cobra round-trip would rewrite the very ids the evaluation aligns on.
    """
    try:
        with open(model_path, encoding="utf-8", errors="replace") as fh:
            head = fh.read(8_000_000)
    except OSError:
        return "cobra"
    if "geneProductAssociation" in head:
        return "cobra"
    if "GENE_ASSOCIATION:" in head:
        return "text"
    return "cobra"


def _or2_find_gff(output_dir, kofam_path):
    """Locate the annotation GFF used only for model-gene/KOfam mapping."""
    candidates = []
    if kofam_path:
        candidates.extend(sorted(glob.glob(os.path.join(os.path.dirname(kofam_path), "*.gff"))))
    candidates.extend(sorted(glob.glob(os.path.join(output_dir, "kofam", "*.gff"))))
    candidates.extend(sorted(glob.glob(os.path.join(output_dir, "*.gff"))))
    return candidates[0] if candidates else ""


def step_gpr_check(output_dir, genome_id, org_name="", taxonomy_tree="",
                   kofam_tsv=None, min_or=2, gene_pos_table=""):
    """Review pure-OR GPRs with the OR2-LLM-v1 protocol.

    Deterministic EC groups are evidence only. The LLM reviews every selected
    branch and returns KEEP/REMOVE/UNCERTAIN. Only REMOVE is applied, a rule is
    never emptied, and benchmark truth is never supplied to the reviewer.
    """
    gapfilled_xml = os.path.join(output_dir, f"{genome_id}-gapfilled.xml")
    if not os.path.exists(gapfilled_xml):
        print("  Skipping GPR review (no gapfilled model).")
        return False

    kofam = kofam_tsv or _find_kofam_tsv(output_dir, genome_id)
    if not os.path.isfile(kofam) or os.path.getsize(kofam) == 0:
        print(f"  Skipping GPR review (no valid KOFam annotation at {kofam}).")
        return False

    gpr_out = os.path.join(output_dir, "gpr_check")
    audit_dir = os.path.join(gpr_out, "audit")
    review_dir = os.path.join(gpr_out, "review")
    os.makedirs(audit_dir, exist_ok=True)
    os.makedirs(review_dir, exist_ok=True)
    checked_xml = os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked.xml")
    pre_sybil = os.path.join(audit_dir, f"{genome_id}-or2-pre-sybil.RDS")
    model_dump = os.path.join(audit_dir, "model_dump.tsv")
    auth_rds = os.path.join(output_dir, f"{genome_id}-gapfilled.RDS")
    if not os.path.exists(auth_rds):
        auth_rds = os.path.join(output_dir, f"{genome_id}-draft.RDS")

    # Build one neutral pre-review representation for reaction names, external
    # notes-style models, and the final model-invariant comparison.
    pre_cmd = [RSCRIPT, f"{CODE_DIR}/agent/agent_gpr_check/export_sybil_model.R",
               "-x", gapfilled_xml, "-o", pre_sybil, "-v"]
    if os.path.exists(auth_rds):
        pre_cmd += ["-m", auth_rds]
    pre_result = run_command(pre_cmd, "OR2: export pre-review sybil model", check=False)
    if (pre_result is None or pre_result.returncode != 0 or
            not os.path.isfile(pre_sybil) or os.path.getsize(pre_sybil) == 0):
        print("  OR2 pre-model export failed.")
        return False

    dump_result = run_command(
        [RSCRIPT, f"{CODE_DIR}/agent/agent_gpr_check/dump_sybil_for_or2.R",
         pre_sybil, model_dump],
        "OR2: dump pre-review GPR table", check=False)
    if dump_result is None or dump_result.returncode != 0 or not os.path.isfile(model_dump):
        print("  OR2 model dump failed.")
        return False

    applier = _gpr_applier(gapfilled_xml)
    gff = _or2_find_gff(output_dir, kofam)
    if applier == "text":
        audit_cmd = [
            PYTHON, f"{CODE_DIR}/agent/agent_gpr_check/audit_or2_gpr_tsv.py",
            "--model_dump", model_dump,
            "--kofam", kofam,
            "--gff", gff,
            "--gene_pos_table", gene_pos_table or "",
            "--out_dir", audit_dir,
            "--min_or", str(min_or),
        ]
    else:
        audit_cmd = [
            PYTHON, f"{CODE_DIR}/agent/agent_gpr_check/gpr_or2_corrected_audit.py",
            "--model", gapfilled_xml,
            "--kofam", kofam,
            "--out_dir", audit_dir,
            "--min_or", str(min_or),
        ]
        if gff:
            audit_cmd += ["--gff", gff]
    audit_result = run_command(audit_cmd, "OR2: build pure-OR evidence table", check=False)
    reactions_tsv = os.path.join(audit_dir, "or2_reactions.tsv")
    genes_tsv = os.path.join(audit_dir, "or2_genes.tsv")
    if (audit_result is None or audit_result.returncode != 0 or
            not os.path.isfile(reactions_tsv) or not os.path.isfile(genes_tsv)):
        print("  OR2 candidate generation failed.")
        return False

    # The dump supplies reaction names. EC clusters in the audit are supporting
    # evidence only; the LLM makes every final branch verdict.
    review_cmd = [
        PYTHON, f"{CODE_DIR}/agent/agent_gpr_check/review_or2_with_llm.py",
        "--audit_dir", audit_dir,
        "--legacy_manifest", model_dump,
        "--kofam", kofam,
        "--out_dir", review_dir,
        "--organism", org_name or genome_id,
        "--cohort", "all",
        "--limit", "0",
        "--workers", str(GPR_WORKERS),
        "--max_retries", str(REVIEW_RETRIES),
        "--llm_config_dir", f"{CODE_DIR}/agent",
    ]
    review_result = run_command(review_cmd, "OR2: LLM review of every pure-OR branch", check=False)
    decisions = os.path.join(review_dir, "llm_decisions.tsv")
    if (review_result is None or review_result.returncode != 0 or
            not os.path.isfile(decisions) or os.path.getsize(decisions) == 0):
        print("  OR2 review failed or produced no decision table.")
        return False

    applier_script = ("apply_or2_llm_text.py" if applier == "text"
                      else "apply_or2_llm_xml.py")
    apply_result = run_command(
        [PYTHON, f"{CODE_DIR}/agent/agent_gpr_check/{applier_script}",
         "--model", gapfilled_xml, "--decisions", decisions, "--out", checked_xml],
        "OR2: apply REMOVE decisions only", check=False)
    if (apply_result is None or apply_result.returncode != 0 or
            not os.path.isfile(checked_xml) or os.path.getsize(checked_xml) == 0):
        print("  OR2 decision application failed.")
        return False

    checked_rds = os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked.RDS")
    write_cmd = [RSCRIPT, f"{CODE_DIR}/agent/agent_gpr_check/write_rds.R",
                 "-x", checked_xml, "-o", checked_rds, "-v"]
    if os.path.exists(auth_rds):
        write_cmd += ["-m", auth_rds]
    write_result = run_command(write_cmd, "OR2: save reviewed RDS", check=False)
    if (write_result is None or write_result.returncode != 0 or
            not os.path.isfile(checked_rds) or os.path.getsize(checked_rds) == 0):
        print("  Saving the OR2-reviewed RDS failed.")
        return False

    _charge_fixer = f"{CODE_DIR}/agent/agent_gpr_check/fix_charges.py"
    if os.path.exists(_charge_fixer):
        run_command([sys.executable, _charge_fixer, checked_xml],
                    "OR2: sanitize SBML charge sentinels", check=False)

    sybil_rds = os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked-sybil.RDS")
    export_cmd = [RSCRIPT, f"{CODE_DIR}/agent/agent_gpr_check/export_sybil_model.R",
                  "-x", checked_xml, "-o", sybil_rds, "-v"]
    maps = [p for p in (checked_rds, auth_rds) if os.path.exists(p)]
    if maps:
        export_cmd += ["-m", ",".join(maps)]
    export_result = run_command(export_cmd, "OR2: export reviewed sybil model", check=False)
    if (export_result is None or export_result.returncode != 0 or
            not os.path.isfile(sybil_rds) or os.path.getsize(sybil_rds) == 0):
        print("  Exporting the OR2-reviewed sybil model failed.")
        return False

    invariant = os.path.join(gpr_out, "invariant.tsv")
    qc_result = run_command(
        [RSCRIPT, f"{CODE_DIR}/agent/agent_gpr_check/check_gpr_only_invariant.R",
         "-b", pre_sybil, "-a", sybil_rds, "-o", invariant],
        "OR2: verify GPR-only model invariants", check=False)
    if qc_result is None or qc_result.returncode != 0 or not os.path.isfile(invariant):
        print("  OR2 invariant check failed to run.")
        return False
    with open(invariant, encoding="utf-8", errors="replace") as fh:
        verdict = [line for line in fh if line.startswith("ALL_PASS")]
    if not verdict or verdict[-1].rstrip().split("\t")[-1] != "YES":
        print("  OR2 model invariant check did not pass; see gpr_check/invariant.tsv")
        return False
    return True


# ==============================================================================
# Main
# ==============================================================================

def step_us_review(output_dir, genome_id, microbe_id, claims_dir, positive_only=False, base="gapseq"):
    """`us-review-expand`: apply positive literature claims to uptake/secretion.

    Applies per-strain claims (<claims_dir>/<microbe_id>.tsv) to the
    GPR-checked model with the canonical ``agent_us_expand`` skill and writes
    <genome_id>-gapfilled-gprchecked-usreviewed.xml (+ .RDS). The stage is a
    no-op (with a warning) when the claims file is absent. Negative claims are
    ignored here; they are handled by ``agent_us_prune``.
    """
    model = os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked.xml")
    if not os.path.isfile(model) or os.path.getsize(model) == 0:
        print(f"  us-review: input model missing: {model}", file=sys.stderr)
        return False
    claims = os.path.join(claims_dir or "", f"{microbe_id}.tsv")
    if not (claims_dir and os.path.isfile(claims) and os.path.getsize(claims) > 0):
        print(f"  us-review: no claims for {microbe_id} ({claims}); stage skipped.")
        return True
    out_xml = os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked-usreviewed.xml")
    out_rds = os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked-usreviewed.RDS")
    eng_dir = os.path.join(output_dir, "us_review_expand")
    os.makedirs(eng_dir, exist_ok=True)
    engine = os.path.join(CODE_DIR, "agent", "agent_us_expand.py")
    cmd = [PYTHON, engine, "--model", model, "--organism", microbe_id,
           "--base", base, "--claims", claims,
           "--out-xml", out_xml, "--out-dir", eng_dir, "--allow-create"]
    exchange_map = os.environ.get("GEMAGENT_US_REVIEW_DIR", "").strip()
    if exchange_map:
        cmd += ["--exchange-map", exchange_map]
    if subprocess.run(cmd).returncode != 0:
        print("  Aborting: us-review expand failed.", file=sys.stderr)
        return False
    if not (os.path.isfile(out_xml) and os.path.getsize(out_xml) > 0):
        print("  Aborting: us-review did not produce a model.", file=sys.stderr)
        return False
    envs = dict(os.environ)
    genv = os.path.join(os.path.dirname(os.path.dirname(PYTHON)), "lib")
    if os.path.isdir(genv):
        envs["LD_LIBRARY_PATH"] = genv + (
            (":" + envs["LD_LIBRARY_PATH"]) if envs.get("LD_LIBRARY_PATH") else "")
    rds_cmd = [RSCRIPT, os.path.join(CODE_DIR, "agent", "agent_gpr_check", "write_rds.R"),
               "-x", out_xml, "-o", out_rds,
               "-m", os.path.join(output_dir, f"{genome_id}-gapfilled-gprchecked.RDS")]
    if subprocess.run(rds_cmd, env=envs).returncode != 0 or not (
            os.path.isfile(out_rds) and os.path.getsize(out_rds) > 0):
        print("  Aborting: us-review RDS export failed.", file=sys.stderr)
        return False
    print(f"  us-review complete: {os.path.basename(out_xml)}")
    return True


def process_genome(input_dir, output_dir, org_name=None,
                   taxonomy_tree=None,
                   reuse_evidence=False, force_rebuild=False, from_draft=False,
                   draft_sbml=None, biomass_template="",
                   gapfill_medium=None, eval_medium=None, kofam_tsv=None,
                   gpr_min_or=2, gene_pos_table=None,
                   gtdb_fallback=True, allow_network=True,
                   biomass_review=True, gapfill_review=True, gpr_check=True,
                   gapfill_review_retry_one_by_one=False,
                   us_review=False, us_claims_dir=None, us_microbe_id=None):
    """Run the complete pipeline for a single genome.

    Front-end = gapseq (find / find-transport / draft / medium / fill).
    Injected LLM reviews = biomass, gapfilling, GPR.

    from_draft=True resumes from an existing draft instead of rebuilding it,
    skipping the gapseq front end (the expensive part); the draft RDS/XML
    must already be present. Use it while iterating on the LLM reviews.
    """

    os.makedirs(output_dir, exist_ok=True)
    os.environ["GEMAGENT_OUTPUT_SUBDIR"] = os.path.basename(output_dir.rstrip("/"))

    try:
        genome_file = find_genome_file(input_dir)
        print(f"\n  Genome file: {genome_file}")
    except FileNotFoundError as e:
        print(f"  Error: {e}", file=sys.stderr)
        return False

    genome_id = extract_genome_id(genome_file)
    print(f"  Genome ID:   {genome_id}")

    # Naming comes from code/taxonomy.py, which owns the whole fallback chain
    # (GTDB r214 table -> GTDB-TK; NCBI RefSeq table -> NCBI API). Two distinct
    # fields on purpose: taxonomy_tree is the GTDB lineage the reviews get as
    # context, org_name is the *NCBI* name that org_prompt_gen matches against
    # MIMEdb -- GTDB renames species MIMEdb does not know, so the two cannot be
    # the same string (see the module docstring).
    tax = taxonomy.resolve(genome_id, genome_file=genome_file,
                           output_dir=output_dir, org_name=org_name or "",
                           taxonomy_tree=taxonomy_tree or "",
                           allow_network=allow_network,
                           gtdb_fallback=gtdb_fallback)
    taxonomy.write_result(output_dir, tax)
    # An explicitly passed name still wins; otherwise use the resolved one.
    # Falling back to genome_id is deliberately last: "bsub" matches nothing in
    # MIMEdb and would silently blank out every review's organism context.
    org_name = org_name or tax.org_name or genome_id
    taxonomy_tree = taxonomy_tree or tax.gtdb_lineage
    print(f"  Organism:    {org_name}")
    if taxonomy_tree:
        print(f"  Taxonomy:    {taxonomy_tree[:80]}...")

    print("\n" + "=" * 60)
    print("  GEMAgent Pipeline (gapseq draft + 3 LLM reviews)")
    print("=" * 60)

    draft_rds_probe = os.path.join(output_dir, f"{genome_id}-draft.RDS")
    draft_xml_probe = os.path.join(output_dir, f"{genome_id}-draft.xml")

    # Cross-tool: seed the draft from an external SBML, then behave exactly like
    # --from_draft (the gapseq front end is skipped; the reviews run on top of
    # the foreign model).
    if draft_sbml:
        if not seed_draft_from_sbml(output_dir, genome_id, draft_sbml):
            return False
        from_draft = True

    if from_draft:
        missing = [p for p in (draft_rds_probe, draft_xml_probe)
                   if not os.path.exists(p)]
        if missing:
            print("  Error: --from_draft needs an existing draft; missing:",
                  file=sys.stderr)
            for _p in missing:
                print(f"    {_p}", file=sys.stderr)
            return False
        print("  Resuming from the existing draft "
              "(skipping gapseq find / find-transport / draft).")
        step_prompts(output_dir, org_name)
    elif not reuse_evidence:
        step_prompts(output_dir, org_name)
        if not step_find_pathways(genome_file, output_dir, genome_id):
            return False
        if not step_predict_transporters(genome_file, output_dir, genome_id):
            return False
        # gapseq draft: native find -> find-transport -> draft.
        if not step_build_draft_gapseq(output_dir, genome_id, genome_file,
                                       force=force_rebuild):
            return False
    else:
        print("  Reusing existing BLAST/pathway/transporter evidence.")
        step_prompts(output_dir, org_name)
        # gapseq draft: native find -> find-transport -> draft.
        if not step_build_draft_gapseq(output_dir, genome_id, genome_file,
                                       force=force_rebuild):
            return False

    # the draft itself records which biomass template gapseq chose; an explicit
    # --biomass_template wins. Foreign models name their biomass reaction
    # differently (ModelSEED: "GramPositiveBiomass_auto_biomass"), so the
    # name-based sniffing in step_read_biomass_template cannot see it and would
    # silently fall back to the Gram_pos default.
    gram_staining = biomass_template or step_read_biomass_template(output_dir, genome_id)
    if biomass_template:
        print(f"  Biomass template (explicit): {biomass_template}")

    # KOFam is evidence for biomass/GPR only. Do not invoke it when both
    # consumers are disabled (e.g. a front-end or gapfill-only run).
    if biomass_review or gpr_check:
        if kofam_tsv:
            if not os.path.isfile(kofam_tsv) or os.path.getsize(kofam_tsv) == 0:
                print(f"  Error: supplied KOFam TSV is missing or empty: {kofam_tsv}", file=sys.stderr)
                return False
            print(f"  Using supplied KOFam annotation: {kofam_tsv}")
        else:
            kofam_tsv = step_annotate_kofam(output_dir, genome_id, genome_file)
            if not kofam_tsv:
                print("  Aborting: enabled agent review requires a valid KOFam annotation.", file=sys.stderr)
                return False

    # growth medium
    #   * the predicted medium is still produced (diagnostics, medium reports),
    #   * but gapfilling runs on _resolve_gapfill_medium(...): the evaluation
    #     medium when one was given via --eval_medium (what the benchmark models
    #     did -- see the policy comment at the top), else the predicted medium.
    if not gapfill_medium:
        if not step_predict_medium(output_dir, genome_id):
            print("  Aborting: growth-medium prediction failed.", file=sys.stderr)
            return False
    gapfill_medium = _resolve_gapfill_medium(output_dir, genome_id, gapfill_medium,
                                             eval_medium)
    if not os.path.isfile(gapfill_medium) or os.path.getsize(gapfill_medium) == 0:
        print(f"  Aborting: gapfill medium is missing or empty: {gapfill_medium}", file=sys.stderr)
        return False
    print(f"  Gapfill medium: {gapfill_medium}")

    # [Review 1] biomass audit -> <id>-draft-biomassadj.RDS
    if biomass_review:
        if not step_biomass_check(
                output_dir, genome_id, gram_staining, org_name=org_name,
                taxonomy_tree=taxonomy_tree or "", force=force_rebuild,
                medium_csv=gapfill_medium, kofam_tsv=kofam_tsv):
            print("  Aborting: biomass review failed.", file=sys.stderr)
            return False
    # step_biomass_check rewrites <id>-draft.RDS *in place* (raw kept as
    # <id>-draft-nobiomass.RDS); resolve the gapfill input via
    # _gapfill_model_file so the name contract with gf.suite.R holds.
    draft_rds = _gapfill_model_file(output_dir, genome_id)
    print(f"  Gapfilling model input: {os.path.basename(draft_rds)}")

    # gapfilling (gapseq gf.suite.R)
    if not step_gapfilling(output_dir, genome_id, gapfill_medium=gapfill_medium,
                           model_rds=draft_rds):
        print("  Aborting: gapfilling failed (no feasible solution).")
        return False

    # Ablation snapshot: keep the pre-review gapfilled model (A2) before the
    # LLM bypass review may re-run gapfilling and overwrite it. Idempotent.
    for _ext in ("RDS", "xml"):
        _src = os.path.join(output_dir, f"{genome_id}-gapfilled.{_ext}")
        _dst = os.path.join(output_dir, f"{genome_id}-gapfilled-prereview.{_ext}")
        if os.path.exists(_src) and not os.path.exists(_dst):
            shutil.copy2(_src, _dst)

    # [Review 2] gapfilling review
    if gapfill_review:
        if not step_review_gapfill_added(
                output_dir, genome_id, org_name=org_name or "",
                taxonomy_tree=taxonomy_tree or "", gapfill_medium=gapfill_medium,
                eval_medium=eval_medium,
                retry_one_by_one=gapfill_review_retry_one_by_one):
            print("  Aborting: gapfill review failed.", file=sys.stderr)
            return False

    # [Review 3] GPR review
    if gpr_check:
        if not step_gpr_check(
                output_dir, genome_id, org_name=org_name or "",
                taxonomy_tree=taxonomy_tree or "", kofam_tsv=kofam_tsv,
                min_or=gpr_min_or, gene_pos_table=gene_pos_table or ""):
            print("  Aborting: GPR review failed.", file=sys.stderr)
            return False

    # [Review 4] us review (literature-driven uptake/secretion curation)
    if us_review:
        if not step_us_review(output_dir, genome_id,
                              us_microbe_id or genome_id,
                              us_claims_dir or "", positive_only=True):
            print("  Aborting: us review failed.", file=sys.stderr)
            return False

    print("\n" + "=" * 60)
    print(f"  Pipeline complete for {genome_id}")
    print("=" * 60)
    return True


def main():
    parser = argparse.ArgumentParser(
        description="GEMAgent Ã¢â‚¬â€ gapseq draft + LLM review pipeline"
    )
    parser.add_argument("--input_dir", type=str, required=True,
                        help="Directory containing genome sequence file")
    parser.add_argument("--output_dir", type=str, default="",
                        help="Output directory (default: input_dir/gemagent_output)")
    parser.add_argument("--reuse_evidence", action="store_true",
                        help="Skip BLAST/transporter if tables already exist")
    parser.add_argument("--force_rebuild", action="store_true",
                        help="Rebuild draft and re-run all reviews")
    parser.add_argument("--from_draft", action="store_true",
                        help="Resume from an existing draft (skips gapseq find/find-transport/draft); errors if it is missing")
    parser.add_argument("--draft_sbml", type=str, default="",
                        help="Seed the draft from an external SBML model (e.g. a "
                             "ModelSEED / CarveMe reconstruction) instead of building "
                             "a gapseq draft; implies --from_draft. The SBML should "
                             "already carry gapseq-style coordinate gene names -- "
                             "see cross_tool/relabel_model.py")
    parser.add_argument("--biomass_template", type=str, default="",
                        choices=["", "Gram_pos", "Gram_neg", "archaea"],
                        help="Explicit biomass template group, overriding the one "
                             "sniffed from the draft. Foreign models name their "
                             "biomass reaction differently, so sniffing silently "
                             "falls back to Gram_pos")
    parser.add_argument("--gene_pos_table", type=str, default="",
                        help="TSV (gene,contig,start,end) bridging foreign, "
                             "non-coordinate gene IDs to genomic coordinates for "
                             "the GPR review. Only needed when the model was NOT "
                             "relabelled with cross_tool/relabel_model.py")
    parser.add_argument("--org_name", type=str, default="",
                        help="Organism name (e.g., 'Escherichia coli'). "
                             "Default: resolved from the genome's assembly "
                             "accession via code/taxonomy.py (NCBI name -- the "
                             "one MIMEdb actually knows)")
    parser.add_argument("--taxonomy_tree", type=str, default="",
                        help="GTDB taxonomy classification string. Default: "
                             "resolved from the accession by code/taxonomy.py")
    parser.add_argument("--no_gtdb_fallback", action="store_true",
                        help="Do not run GTDB-TK when the offline GTDB r214 "
                             "table has no entry for this genome (that run "
                             "costs 10-30 min)")
    parser.add_argument("--no_network", action="store_true",
                        help="Offline mode: skip the NCBI Datasets API lookup "
                             "when the local RefSeq table has no entry")
    parser.add_argument("--gapfill_medium", type=str, default="",
                        help="Explicit medium CSV for gapfilling (overrides "
                             "--gapfill_medium_mode)")
    parser.add_argument("--eval_medium", type=str, default="",
                        help="The medium CSV this model will be evaluated on. "
                             "Gapfilling then runs on it (so the model is "
                             "guaranteed to grow where its F1 is measured). "
                             "This is what the benchmark models did.")
    parser.add_argument("--gapfill_medium_mode", type=str, default="",
                        choices=["eval", "complete", "predicted", ""],
                        help="Gapfill medium policy when --gapfill_medium is not "
                             "given: 'eval' = --eval_medium (default when that is "
                             "given), 'predicted' = the organism's own predicted "
                             "medium, 'complete' = gapseq's ALLmed.csv "
                             "(experimental; NOT what the benchmark used)")
    parser.add_argument("--kofam_tsv", type=str, default="",
                        help="KOFamScan annotation TSV "
                             "(default: <output_dir>/kofam/<genome_id>_kofam.tsv)")
    parser.add_argument("--gpr_min_or", type=int, default=2,
                        help="Minimum OR-branch count to review a GPR (default 2)")
    parser.add_argument("--no_biomass_review", action="store_true",
                        help="Skip the LLM biomass review")
    parser.add_argument("--gapfill_review_retry_one_by_one", action="store_true",
                        help="Force the gapfill bypass review to retry exclusions "
                             "one-by-one instead of reverting the whole batch")
    parser.add_argument("--no_gapfill_review", action="store_true",
                        help="Skip the post-gapfill LLM bypass review")
    parser.add_argument("--no_gpr_check", action="store_true",
                        help="Skip the post-gapfill LLM GPR review")
    parser.add_argument("--us_claims_dir", type=str, default="",
                        help="Directory with per-strain us-review claims "
                             "(<MicrobeID>.tsv); enables the us-review stage")
    parser.add_argument("--us_microbe_id", type=str, default="",
                        help="Panel MicrobeID for the us-review claims lookup "
                             "(default: genome_id)")
    parser.add_argument("--no_us_review", action="store_true",
                        help="Skip the post-GPR us review")
    args = parser.parse_args()

    if args.gapfill_medium_mode:
        global GAPFILL_MEDIUM_MODE
        GAPFILL_MEDIUM_MODE = args.gapfill_medium_mode
        print(f"  Gapfill medium mode: {GAPFILL_MEDIUM_MODE}")

    if not os.path.isdir(args.input_dir):
        print(f"Error: input directory not found: {args.input_dir}", file=sys.stderr)
        sys.exit(1)

    output_dir = args.output_dir or os.path.join(args.input_dir, OUTPUT_SUBDIR)

    success = process_genome(
        input_dir=args.input_dir,
        output_dir=output_dir,
        org_name=args.org_name if args.org_name else None,
        taxonomy_tree=args.taxonomy_tree if args.taxonomy_tree else None,
        reuse_evidence=args.reuse_evidence,
        force_rebuild=args.force_rebuild,
        from_draft=args.from_draft,
        draft_sbml=args.draft_sbml if args.draft_sbml else None,
        biomass_template=args.biomass_template,
        gapfill_medium=args.gapfill_medium if args.gapfill_medium else None,
        eval_medium=args.eval_medium if args.eval_medium else None,
        kofam_tsv=args.kofam_tsv if args.kofam_tsv else None,
        gpr_min_or=args.gpr_min_or,
        gene_pos_table=args.gene_pos_table if args.gene_pos_table else None,
        gtdb_fallback=not args.no_gtdb_fallback,
        allow_network=not args.no_network,
        biomass_review=not args.no_biomass_review,
        gapfill_review=not args.no_gapfill_review,
        gapfill_review_retry_one_by_one=args.gapfill_review_retry_one_by_one,
        gpr_check=not args.no_gpr_check,
        us_review=(bool(args.us_claims_dir) and not args.no_us_review),
        us_claims_dir=args.us_claims_dir if args.us_claims_dir else None,
        us_microbe_id=args.us_microbe_id if args.us_microbe_id else None,
    )
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
