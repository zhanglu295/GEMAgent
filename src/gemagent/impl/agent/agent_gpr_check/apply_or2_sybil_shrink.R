#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(sybil)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("usage: apply_or2_sybil_shrink.R <base.RDS> <actions.tsv> <out.RDS> <report.tsv>")
}

base_path <- args[[1]]
actions_path <- args[[2]]
out_path <- args[[3]]
report_path <- args[[4]]
sybil::SYBIL_SETTINGS("SOLVER", "glpkAPI")

model <- readRDS(base_path)
actions <- read.delim(actions_path, stringsAsFactors = FALSE, quote = "")
if (!all(c("reaction_id", "gpr_after") %in% names(actions))) {
  stop("actions file must contain reaction_id and gpr_after")
}

growth_before <- tryCatch(
  as.numeric(optimizeProb(model)@lp_obj),
  error = function(e) NA_real_
)
genes_before <- as.character(model@allGenes)
applied <- 0L
missing_reactions <- character()

for (i in seq_len(nrow(actions))) {
  rid <- actions$reaction_id[[i]]
  ri <- match(rid, model@react_id)
  if (is.na(ri)) {
    missing_reactions <- c(missing_reactions, rid)
    next
  }
  rule0 <- actions$gpr_after[[i]]
  if (is.na(rule0) || !nzchar(rule0)) {
    stop(sprintf("empty gpr_after for %s", rid))
  }
  tokens <- regmatches(rule0, gregexpr("[A-Za-z0-9_.:+-]+", rule0))[[1]]
  tokens <- unique(tokens[!grepl("^(and|or)$", tokens, ignore.case = TRUE)])
  missing_genes <- setdiff(tokens, model@allGenes)
  if (length(missing_genes)) {
    stop(sprintf("unmapped genes for %s: %s", rid, paste(missing_genes, collapse = ",")))
  }

  rule <- rule0
  for (gene in tokens[order(nchar(tokens), decreasing = TRUE)]) {
    rule <- gsub(gene, paste0("x[", match(gene, tokens), "]"), rule, fixed = TRUE)
  }
  rule <- gsub(" and ", " & ", rule, fixed = TRUE)
  rule <- gsub(" or ", " | ", rule, fixed = TRUE)

  model@genes[[ri]] <- tokens
  model@gprRules[ri] <- rule
  model@gpr[ri] <- gsub(" and ", " & ", rule0, fixed = TRUE) |>
    gsub(pattern = " or ", replacement = " | ", fixed = TRUE)
  applied <- applied + 1L
}

used_after <- unique(unlist(model@genes, use.names = FALSE))
genes_after <- genes_before[genes_before %in% used_after]
orphaned_removed <- setdiff(genes_before, genes_after)

model@allGenes <- genes_after
gene_index <- setNames(seq_along(genes_after), genes_after)
mat <- Matrix::Matrix(
  0,
  nrow = length(model@react_id),
  ncol = length(genes_after),
  sparse = TRUE
)
for (ri in seq_along(model@genes)) {
  genes <- unique(model@genes[[ri]])
  genes <- genes[genes %in% genes_after]
  model@genes[ri] <- list(genes)
  if (length(genes)) {
    mat[ri, unname(gene_index[genes])] <- 1
  }
}
model@rxnGeneMat <- mat

growth_after <- tryCatch(
  as.numeric(optimizeProb(model)@lp_obj),
  error = function(e) NA_real_
)
if (is.na(growth_before) || is.na(growth_after) ||
    growth_after < growth_before - max(1e-9, growth_before * 1e-6)) {
  stop(sprintf("growth invariant failed: %.12f -> %.12f", growth_before, growth_after))
}

saveRDS(model, out_path)
report <- data.frame(
  applied = applied,
  reactions = length(model@react_id),
  metabolites = length(model@met_id),
  genes_before = length(genes_before),
  genes_after = length(genes_after),
  orphaned_genes_removed = length(orphaned_removed),
  growth_before = growth_before,
  growth_after = growth_after,
  missing_reactions = paste(missing_reactions, collapse = ";"),
  model = out_path,
  stringsAsFactors = FALSE
)
write.table(report, report_path, sep = "\t", quote = FALSE, row.names = FALSE)
print(report)
cat("APPLY_OR2_SYBIL_SHRINK_DONE\n")
