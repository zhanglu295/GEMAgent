#!/usr/bin/env Rscript
# Export metadata for reactions added by official gapseq fill.
# This is a reporting helper, not a fork of gf.suite.R.
suppressPackageStartupMessages({library(data.table); library(jsonlite)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) stop("usage: gapfill_evidence.R draft_rds gapfilled_rds weights_rds genes_rds seed_db_tsv out_json")
draft_rds <- args[[1]]; filled_rds <- args[[2]]; weights_rds <- args[[3]]
genes_rds <- args[[4]]; seed_tsv <- args[[5]]; out_json <- args[[6]]
mod0 <- readRDS(draft_rds); mod1 <- readRDS(filled_rds)
ids0 <- gsub("_c0$", "", mod0@react_id)
ids1 <- gsub("_c0$", "", mod1@react_id)
added <- sort(unique(setdiff(ids1, ids0)))
added <- added[!grepl("^EX_", added)]
w <- as.data.table(readRDS(weights_rds)); g <- as.data.table(readRDS(genes_rds))
seeddb <- if (file.exists(seed_tsv)) fread(seed_tsv) else data.table()
rows <- lapply(added, function(sid) {
  db <- if ("id" %in% names(seeddb)) seeddb[id == sid] else seeddb[0]
  wr <- if ("seed" %in% names(w)) w[seed == sid] else w[0]
  gr <- if ("seed" %in% names(g)) g[seed == sid] else g[0]
  if ("rm" %in% names(gr)) gr <- gr[is.na(rm) | rm == FALSE]
  genes <- if ("gene" %in% names(gr)) unique(as.character(gr$gene)) else character()
  genes <- genes[!is.na(genes) & nzchar(genes)]
  list(
    seed = sid,
    name = if (nrow(db) && "name" %in% names(db)) as.character(db$name[1]) else "",
    equation = if (nrow(db) && "definition" %in% names(db)) as.character(db$definition[1]) else "",
    ec = if (nrow(db) && "ec_numbers" %in% names(db)) as.character(db$ec_numbers[1]) else "",
    pathway = if (nrow(wr) && "pathway" %in% names(wr)) as.character(wr$pathway[1]) else "",
    pathway_status = if (nrow(wr) && "pathway.status" %in% names(wr)) as.character(wr$pathway.status[1]) else "",
    src = if (nrow(wr) && "src" %in% names(wr)) as.character(wr$src[1]) else "",
    pident = if (nrow(wr) && "pident" %in% names(wr)) suppressWarnings(max(wr$pident, na.rm = TRUE)) else NA_real_,
    bitscore = if (nrow(wr) && "bitscore" %in% names(wr)) suppressWarnings(max(wr$bitscore, na.rm = TRUE)) else NA_real_,
    genes = paste(genes, collapse = "|"),
    geneless = length(genes) == 0,
    spontaneous = if (nrow(wr) && "spont" %in% names(wr)) isTRUE(as.logical(wr$spont[1])) else FALSE
  )
})
write_json(if (length(rows)) rows else list(), out_json, auto_unbox = TRUE, pretty = TRUE)