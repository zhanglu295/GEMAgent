#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(sybil))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: dump_sybil_for_or2.R <model.RDS> <out.tsv>")
model <- readRDS(args[[1]])
out <- data.frame(
  reaction_id = as.character(model@react_id),
  name = as.character(model@react_name),
  gpr = as.character(model@gpr),
  stringsAsFactors = FALSE
)
write.table(out, args[[2]], sep = "\t", quote = FALSE, row.names = FALSE)
