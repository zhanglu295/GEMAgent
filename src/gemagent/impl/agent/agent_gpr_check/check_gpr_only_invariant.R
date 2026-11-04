#!/usr/bin/env Rscript
# Validate that GPR review changed only gene-protein-reaction annotations.
#
# This is a contract test for the evaluated sybil models, not just for the
# intermediate SBML files. Two models must have identical stoichiometry,
# bounds, objective coefficients and FBA growth. The post-review model may
# drop genes, but it may not introduce genes, reactions or metabolites.

suppressPackageStartupMessages({
  library(getopt)
  library(Matrix)
  library(sybil)
})

spec <- matrix(c(
  "before", "b", 1, "character", "Pre-review sybil model RDS",
  "after",  "a", 1, "character", "Post-review sybil model RDS",
  "out",    "o", 1, "character", "Two-column TSV audit output",
  "tol",    "t", 1, "double",    "Numeric tolerance for invariant checks",
  "verbose","v", 0, "logical",   "Print audit output to stderr",
  "help",   "h", 0, "logical",   "Show this help"
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)
if (!is.null(opt$help) || is.null(opt$before) || is.null(opt$after) ||
    is.null(opt$out)) {
  cat(getopt(spec, usage = TRUE))
  quit(status = 1)
}

if (is.null(opt$tol)) opt$tol <- 1e-9
if (is.null(opt$verbose)) opt$verbose <- FALSE
if (!is.finite(opt$tol) || opt$tol < 0) stop("--tol must be finite and >= 0")

read_model <- function(path, label) {
  if (!file.exists(path)) stop(sprintf("%s model does not exist: %s", label, path))
  model <- readRDS(path)
  if (!methods::is(model, "modelorg")) {
    stop(sprintf("%s object is not a sybil modelorg: %s", label, path))
  }
  model
}

before <- read_model(opt$before, "before")
after <- read_model(opt$after, "after")

fmt_num <- function(x) {
  if (length(x) == 0) return("NA")
  if (is.na(x)) return("NA")
  formatC(as.numeric(x), digits = 10, format = "fg", flag = "#")
}

safe_equal <- function(x, y) {
  isTRUE(all.equal(x, y, tolerance = opt$tol, check.attributes = FALSE))
}

same_set <- function(x, y) {
  setequal(as.character(x), as.character(y)) &&
    length(x) == length(y) &&
    !anyDuplicated(x) &&
    !anyDuplicated(y)
}

fba_growth <- function(model) {
  if (!any(model@obj_coef != 0)) return(NA_real_)
  tryCatch(
    as.numeric(sybil::optimizeProb(model)@lp_obj),
    error = function(e) NA_real_
  )
}

normalize_gene <- function(x) {
  # ModelSEED's SBML writer removes its display-only "G_" namespace prefix on
  # round-trip. Normalizing that one representation detail keeps the audit
  # focused on whether genes were incorrectly introduced or removed.
  sub("^G_", "", as.character(x))
}

n_gpr <- function(model) {
  if (!("gpr" %in% methods::slotNames(model))) return(NA_integer_)
  sum(nzchar(trimws(as.character(model@gpr))))
}

matrix_diff_nnz <- function(x, y, tol) {
  if (!identical(dim(x), dim(y))) return(Inf)
  delta <- x - y
  Matrix::nnzero(abs(delta) > tol)
}

react_before <- as.character(before@react_id)
react_after <- as.character(after@react_id)
met_before <- as.character(before@met_id)
met_after <- as.character(after@met_id)
gene_before <- as.character(before@allGenes)
gene_after <- as.character(after@allGenes)

growth_before <- fba_growth(before)
growth_after <- fba_growth(after)

checks <- list()
checks$check_react_id_set_equal <- same_set(react_before, react_after)
checks$check_met_id_set_equal <- same_set(met_before, met_after)
checks$check_S_identical <- identical(dim(before@S), dim(after@S)) &&
  matrix_diff_nnz(before@S, after@S, opt$tol) == 0L
checks$check_nnz_equal <- Matrix::nnzero(before@S) == Matrix::nnzero(after@S)
checks$check_lowbnd_identical <- safe_equal(before@lowbnd, after@lowbnd)
checks$check_uppbnd_identical <- safe_equal(before@uppbnd, after@uppbnd)
checks$check_obj_coef_identical <- safe_equal(before@obj_coef, after@obj_coef)
checks$check_growth_equal <- !is.na(growth_before) &&
  !is.na(growth_after) &&
  abs(growth_before - growth_after) <= max(opt$tol, opt$tol * abs(growth_before))
checks$check_genes_monotone <- all(normalize_gene(gene_after) %in% normalize_gene(gene_before)) &&
  !anyDuplicated(normalize_gene(gene_after))
checks$check_gpr_nonempty <- is.finite(n_gpr(after)) && n_gpr(after) > 0L

all_pass <- all(vapply(checks, isTRUE, logical(1)))

rows <- list(
  c("before", opt$before),
  c("after", opt$after),
  c("n_react_before", length(react_before)),
  c("n_react_after", length(react_after)),
  c("n_met_before", length(met_before)),
  c("n_met_after", length(met_after)),
  c("n_genes_before", length(gene_before)),
  c("n_genes_after", length(gene_after)),
  c("growth_before", fmt_num(growth_before)),
  c("growth_after", fmt_num(growth_after)),
  c("nnz_before", Matrix::nnzero(before@S)),
  c("nnz_after", Matrix::nnzero(after@S)),
  c("n_gpr_before", n_gpr(before)),
  c("n_gpr_after", n_gpr(after))
)

for (name in names(checks)) {
  rows[[length(rows) + 1L]] <- c(name, if (checks[[name]]) "PASS" else "FAIL")
}
rows[[length(rows) + 1L]] <- c("ALL_PASS", if (all_pass) "YES" else "NO")

audit <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
names(audit) <- c("metric", "value")
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
write.table(
  audit, file = opt$out, sep = "\t", quote = FALSE,
  row.names = FALSE, col.names = TRUE
)

if (opt$verbose) {
  print(audit, row.names = FALSE)
}
quit(status = if (all_pass) 0 else 2)
