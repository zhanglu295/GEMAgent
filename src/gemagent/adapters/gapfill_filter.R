#!/usr/bin/env Rscript
# Filter gapseq candidate tables by a review exclusion list.
# This is not a gapfilling implementation; it only removes selected seeds
# before calling the official `gapseq fill` CLI.
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("usage: gapfill_filter.R weights_in genes_in exclude_tsv weights_out genes_out")
weights_in <- args[[1]]; genes_in <- args[[2]]; exclude <- args[[3]]
weights_out <- args[[4]]; genes_out <- args[[5]]
w <- readRDS(weights_in)
g <- readRDS(genes_in)
ex <- if (file.exists(exclude)) fread(exclude) else data.table(seed = character())
seeds <- if ("seed" %in% names(ex)) as.character(ex$seed) else character()
seeds <- seeds[!is.na(seeds) & nzchar(seeds)]
if (length(seeds) && "seed" %in% names(w)) {
  if (data.table::is.data.table(w)) w <- w[!(seed %in% seeds)] else w <- w[!(w$seed %in% seeds), , drop = FALSE]
}
if (length(seeds) && "seed" %in% names(g)) {
  if (data.table::is.data.table(g)) g <- g[!(seed %in% seeds)] else g <- g[!(g$seed %in% seeds), , drop = FALSE]
}
saveRDS(w, weights_out)
saveRDS(g, genes_out)