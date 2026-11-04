#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Save a GEMAgent SBML model as an RDS ModelOrg object.
#
# Why this exists
# ---------------
# gapseq persists every reconstruction stage as an RDS object
# (code/gapfilling/gf.suite.R:941-944):
#
#     out.rds <- paste0(output.dir,"/",out.id,"-gapfilled.RDS")
#     saveRDS(mod.out, file = out.rds)
#
# The GPR-review stage, however, is implemented in Python
# (agent_gpr_check/apply_or2_llm_xml.py) and only wrote SBML, so the
# pipeline ended without a post-GPR RDS. This script closes that gap by
# reading the reviewed SBML back into a `cobrar` ModelOrg and saving it
# with exactly the same idiom gapseq uses.
#
# Gene names
# ----------
# This script reads the model back from SBML, which rewrites "." and ":" in
# gene names to "_" and therefore destroys the contig/coordinate structure.
# Pass the pre-GPR RDS via -m and the names are restored (see
# restore_gene_names.R for the full explanation). Without -m they stay in
# their SBML-sanitized form.
#
# Usage
# -----
#   Rscript write_rds.R -x <model.xml> -o <model.RDS> [-m <auth.RDS>] [-v]
# ---------------------------------------------------------------------------
library(getopt)

spec <- matrix(c(
  'xml',     'x', 1, "character", "SBML model to convert",
  'out',     'o', 1, "character", "Output RDS path",
  'map',     'm', 1, "character", "Authoritative pre-GPR RDS with intact (gapseq-style) gene names; used to undo the SBML id sanitization",
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

# ---- Restore gene names destroyed by the SBML id sanitization --------------
if (!is.null(opt$map)) {
  if (!file.exists(opt$map)) stop(paste("Map RDS not found:", opt$map))
  lookup <- gene_name_lookup(readRDS(opt$map)@allGenes)
  attr(lookup, "source") <- opt$map
  mod <- restore_gene_names(mod, lookup, verbose = opt$verbose)
}

# Restore the biomass objective in case an upstream stage left the temporary
# gapfilling target (cpd11416) in place -- same guard as gf.suite.R:932-939.
if ("bio1" %in% mod@react_id) {
  mod@obj_coef <- rep(0, cobrar::react_num(mod))
  mod@obj_coef[which(mod@react_id == "bio1")] <- 1
} else {
  warning("bio1 not found in model; SBML objective may be incorrect.")
}

saveRDS(mod, file = opt$out)

if (opt$verbose) {
  cat(sprintf("  write_rds: react=%d met=%d genes=%d EX=%d -> %s\n",
              cobrar::react_num(mod), cobrar::met_num(mod),
              length(mod@allGenes),
              length(grep("^EX_", mod@react_id)), opt$out))
}