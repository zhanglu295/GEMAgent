#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Export a GEMAgent SBML model as a *sybil* `modelorg` RDS.
#
# Why this exists
# ---------------
# gapseq 2.1 builds its models with `cobrar`, so every GEMAgent-v2 stage
# (draft / gapfilled / gprchecked) saves a `cobrar::ModelOrg` RDS. The official
# gapseqEval (compare_GE_recons.R / predGeneEssentiality.GS.R) is written
# against *sybil*, and `optimizeProb` does not dispatch on `ModelOrg`:
#
#     unable to find an inherited method for function 'optimizeProb'
#     for signature 'object = "ModelOrg"'
#
# The only lossless bridge is SBML -- `sybilSBML::readSBMLmod` reconstructs a
# proper sybil `modelorg` (S matrix, bounds, objective, GPRs all verified to
# give the identical FBA objective). SBML cannot carry "." or ":" in ids
# though, so gene names come back sanitized; -m restores them from an RDS that
# still has them intact (see restore_gene_names.R).
#
# Usage
# -----
#   Rscript export_sybil_model.R -x <model.xml> -o <model-sybil.RDS> \
#                                [-m <auth.RDS> | -a <auth.RDS>[,<auth.RDS>...]] [-v]
# ---------------------------------------------------------------------------
library(getopt)

spec <- matrix(c(
  'xml',     'x', 1, "character", "GEMAgent SBML model to convert",
  'out',     'o', 1, "character", "Output sybil modelorg RDS path",
  'map',     'm', 1, "character", "Authoritative RDS (or comma-separated list) with intact gapseq-style gene names",
  'verbose', 'v', 0, "logical",   "Verbose output. Default: FALSE",
  'help',    'h', 0, "logical",   "help"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

if (!is.null(opt$help) || is.null(opt$xml) || is.null(opt$out)) {
  cat(getopt(spec, usage = TRUE))
  q(status = 1)
}
if (is.null(opt$verbose)) opt$verbose <- FALSE

suppressMessages({ library(sybil); library(sybilSBML) })

# 与 code/gapfilling/*.R 一致的 script.dir 惯用法
initial.options <- commandArgs(trailingOnly = FALSE)
script.name <- sub("--file=", "", initial.options[grep("--file=", initial.options)])
script.dir  <- dirname(script.name)
source(file.path(script.dir, "restore_gene_names.R"))

if (!file.exists(opt$xml)) stop(paste("Model file not found:", opt$xml))

mod <- sybilSBML::readSBMLmod(opt$xml)   # -> sybil modelorg

# ---- Restore gene names mangled by the SBML id sanitization ----------------
if (!is.null(opt$map)) {
  n.before <- length(mod@allGenes)
  lookup <- gene_name_lookup_from_rds(strsplit(opt$map, ",")[[1]])
  if (is.null(lookup)) {
    warning("none of the -m candidate RDS files exist; gene names stay sanitized")
  } else {
    if (length(lookup) < n.before)
      warning(sprintf(paste0("gene-name authority covers only %d names but the ",
                             "model has %d genes; the remainder stay sanitized"),
                      length(lookup), n.before))
    mod <- restore_gene_names(mod, lookup, verbose = opt$verbose)
  }
}

# FBA sanity check: a model that cannot grow is almost always broken wiring.
if ("bio1" %in% mod@react_id) {
  mod@obj_coef <- rep(0, sybil::react_num(mod))
  mod@obj_coef[which(mod@react_id == "bio1")] <- 1
  grow <- tryCatch(as.numeric(sybil::optimizeProb(mod)@lp_obj), error = function(e) NA_real_)
  if (opt$verbose) cat(sprintf("  FBA (unconstrained) objective = %.6f\n", grow))
  if (!is.na(grow) && grow < 1e-6)
    warning("biomass flux is ~0 in the exported model; check the upstream stages")
}

saveRDS(mod, file = opt$out)

if (opt$verbose) {
  n.parseable <- sum(grepl("_[0-9]+:[0-9]+$", mod@allGenes))
  cat(sprintf("  export_sybil_model: class=%s react=%d met=%d genes=%d EX=%d\n",
              paste(class(mod), collapse = "/"), sybil::react_num(mod),
              sybil::met_num(mod), length(mod@allGenes),
              length(grep("^EX_", mod@react_id))))
  cat(sprintf("    coordinate-parseable gene names: %d/%d\n",
              n.parseable, length(mod@allGenes)))
  cat(sprintf("    -> %s\n", opt$out))
}