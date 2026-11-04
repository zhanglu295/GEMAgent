#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Seed a GEMAgent draft RDS (`cobrar` ModelOrg) from an external SBML model.
#
# Why this exists
# ---------------
# GEMAgent's gapfilling (code/gapfilling/gf.suite.R) and the biomass-review
# apply step both consume a `cobrar` ModelOrg RDS (<id>-draft.RDS), not SBML.
# To graft the LLM reviews onto a model built by another tool (ModelSEED,
# CarveMe, ...) we therefore have to materialise that RDS from the foreign
# SBML.
#
# `write_rds.R` cannot be reused verbatim: it reads with cobrar, which mangles
# characters that are illegal in an SBML SId. Three different escapes show up
# across the readers/writers in play:
#
#     "_"  substitution        NZ_CP010052_1_410_1751      (cobrar writer)
#     "__<ord>__"              NZ_CP010052__46__1_410__58__1751 (python-cobra)
#     "gp_...__NAME__"         gp_NZ_CP__ZERO__...__DOT__...
#                              __COLON__                (cobrar reader, 实测)
#
# The third is what cobrar::readSBMLmod actually produces for the external
# drafts, and it is the reason the coordinate structure used to be lost
# entirely (0/855 parseable). decode_sbml_id_any() in restore_gene_names.R
# undoes whichever one is in play.
#
# Every downstream consumer parses the contig/coordinate structure back out of
# the gene name -- the GPR review's KOFam bridge (gpr_check_prompt.py) and the
# official gapseqEval `predGeneEssentiality.GS`. Left encoded, they see a
# degenerate 1-bp range for every gene: KO coverage collapses to 0 and every
# knockout becomes a no-op. So the names are decoded back here, before saving.
#
# The input SBML is expected to carry gapseq-style coordinate names already
# (see cross_tool/relabel_model.py); this script only undoes the SBML-id escape.
#
# Usage
# -----
#   Rscript seed_foreign_draft.R -x <model.xml> -o <draft.RDS> [-v]
# ---------------------------------------------------------------------------
library(getopt)

spec <- matrix(c(
  'xml',     'x', 1, "character", "External SBML model to seed the draft from",
  'out',     'o', 1, "character", "Output cobrar ModelOrg RDS path",
  'verbose', 'v', 0, "logical",   "Verbose output. Default: FALSE",
  'help',    'h', 0, "logical",   "help"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$xml) || is.null(opt$out)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}
if (is.null(opt$verbose)) opt$verbose <- FALSE

suppressMessages(library(cobrar))

# 与 code/gapfilling/*.R 一致的 script.dir 惯用法
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.dir  <- dirname(script.name)
source(file.path(script.dir, "restore_gene_names.R"))

if (!file.exists(opt$xml)) stop(paste("Model file not found:", opt$xml))

mod <- cobrar::readSBMLmod(opt$xml)

# ---- Undo the SBML-id escape on gene names ---------------------------------
# cobrar::readSBMLmod 对同一个基因名有两套（互不相通的）转义：早期只用
# "__<ord>__"（python-cobra 的写法），现在实测外源 SBML 走的是
# "gp_...__DOT__/__COLON__/__ZERO__..."（见 restore_gene_names.R 的
# COBRAR_CHAR_TOKEN 注释）。decode_sbml_id_any() 依次套两种解码，
# 任何一套命中都能把坐标还原，避免 draft 里全是 "gp_..." 导致
# GPR review 的 KOFam 桥与 GE 评测同时失效。
n.before <- length(mod@allGenes)
raw.names <- mod@allGenes
mod@allGenes <- decode_sbml_id_any(mod@allGenes)
mod@genes <- lapply(mod@genes, function(v) {
  if (length(v) == 0) return(v)
  decode_sbml_id_any(v)
})

n.decoded <- sum(mod@allGenes != raw.names)
n.parseable <- sum(grepl("_[0-9]+:[0-9]+$", mod@allGenes))
if (opt$verbose) {
  cat(sprintf("  draft seed: react=%d met=%d genes=%d EX=%d\n",
              cobrar::react_num(mod), cobrar::met_num(mod),
              length(mod@allGenes), length(grep("^EX_", mod@react_id))))
  cat(sprintf("    gene names decoded from SBML-id escape: %d/%d\n",
              n.decoded, n.before))
  cat(sprintf("    coordinate-parseable gene names: %d/%d\n",
              n.parseable, length(mod@allGenes)))
}
# 外源模型若本来就用非坐标基因名（CarveMe 的 G_NP_391348_1），0 命中是正常的，
# 由 --gene_pos_table 桥接；但若一个都没解出来，几乎必然是转义没还原，
# GPR review 与 GE 评测都会静默失效，这里必须吵出来。
if (n.decoded == 0 && n.before > 0 && all(grepl("^(gp_|G_|__)", mod@allGenes))) {
  warning(sprintf(paste0("seed_foreign_draft.R: %d 个基因名一个都没解出来，",
                         "SBML id 转义可能未被识别（例: %s）"),
                  n.before, paste(head(mod@allGenes, 2), collapse = ", ")),
          call. = FALSE)
}

# ---- Speak cobrar/gapseq's metabolite-id dialect ---------------------------
# Every model gapseq and GEMAgent build names its metabolites `<cpd>[<comp>]`,
# e.g. `cpd00001[c0]`. SBML cannot carry "[" or "]" in an SId, so the writer
# emits `cpd00001_c0` and a reader hands that straight back. A draft seeded from
# a foreign SBML therefore ends up in the underscore dialect while every
# model-side consumer in the pipeline hard-codes the bracket dialect:
#
#   code/draft_construction/replace_biomass.R   addReact(met = dt_bm$id)
#   code/gapfilling/add_missing_exRxns.R        add_met_sink()
#   code/gapfilling/add_missing_exRxns.R        add_missing_exchanges()
#
# Mismatched, the biomass swap silently *adds* ~96 duplicate species instead of
# reusing the existing ones: bio1 then demands metabolites nothing produces and
# its product cpd11416 has no consumer, so growth collapses to 0 and gf.suite.R
# aborts with "cannot produce enough target even when all candidate reactions
# are added" (实测 ModelSEED: met 1490 -> 1586, growth 0.285277 -> 0, gapfill
# 全段 Final growth rate: 0).
#
# The rewrite is a pure relabelling -- S, bounds, GPRs and reaction ids are
# untouched -- and cobrar's SBML writer turns the brackets back into `_c0` on
# the way out, so the evaluation side (sybil, media matching) is unaffected.
normalize_met_ids_to_bracket <- function(mod, verbose = FALSE) {
  if (length(mod@met_id) == 0) return(mod)
  if (any(grepl("[", mod@met_id, fixed = TRUE))) {
    if (verbose) cat("    met ids already in bracket dialect; left as is\n")
    return(mod)
  }
  # cobrar stores the compartment as a character vector parallel to met_id, so
  # the suffix to translate is exactly paste0("_", met_comp) -- never a blind
  # sub() that could eat an underscore inside the compound name.
  comp   <- as.character(mod@met_comp)
  suffix <- paste0("_", comp)
  hit    <- substr(mod@met_id, nchar(mod@met_id) - nchar(suffix) + 1,
                   nchar(mod@met_id)) == suffix
  newid  <- mod@met_id
  newid[hit] <- paste0(substr(mod@met_id[hit], 1,
                              nchar(mod@met_id[hit]) - nchar(suffix[hit])),
                       "[", comp[hit], "]")
  if (anyDuplicated(newid))
    stop("normalize_met_ids_to_bracket(): the relabelling collides with an ",
         "existing metabolite id -- refusing to continue.")
  if (sum(hit) == 0)
    stop(paste0("normalize_met_ids_to_bracket(): no metabolite id ends in ",
                paste0("_", paste(unique(comp), collapse = "/")),
                " -- the foreign model's compartment convention is not ",
                "understood. Downstream biomass/gapfill steps hard-code the ",
                "bracket dialect, so continuing would silently produce a ",
                "model that cannot grow."))
  if (is.data.frame(mod@met_attr) && "id" %in% colnames(mod@met_attr))
    mod@met_attr$id <- newid
  mod@met_id <- newid
  if (verbose)
    cat(sprintf("    met ids -> bracket dialect: %d/%d (compartments: %s)\n",
                sum(hit), length(newid), paste(unique(comp), collapse = ",")))
  mod
}

mod <- normalize_met_ids_to_bracket(mod, verbose = opt$verbose)
if (opt$verbose && any(grepl("[", mod@met_id, fixed = TRUE))) {
  ex <- grep("^EX_", mod@react_id, value = TRUE)
  cat(sprintf("    EX_ reactions: %d (e.g. %s)\n", length(ex),
              paste(head(ex, 3), collapse = ", ")))
}

# Restore the biomass objective in case an upstream stage left the temporary
# gapfilling target (cpd11416) in place -- same guard as write_rds.R / gf.suite.R.
if ("bio1" %in% mod@react_id) {
  mod@obj_coef <- rep(0, cobrar::react_num(mod))
  mod@obj_coef[which(mod@react_id == "bio1")] <- 1
} else {
  warning("bio1 not found in model; SBML objective may be incorrect.")
}

saveRDS(mod, file = opt$out)

if (opt$verbose) cat(sprintf("    -> %s\n", opt$out))
