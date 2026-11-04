# ---------------------------------------------------------------------------
# restore_gene_names.R
# ----------------------------------------------------------------------------
# Undo the gene-name mangling that an SBML round-trip inflicts.
#
# SBML `SId` must match [A-Za-z_][A-Za-z0-9_]*, so the SBML writer rewrites
# every gene name by replacing each illegal character with "_":
#
#     gapseq style    NZ_CP010052.1_4168217:4169140
#     after SBML      NZ_CP010052_1_4168217_4169140     ( "." and ":" -> "_" )
#
# The dot and the colon are exactly the two characters that carry the
# contig / coordinate structure. Any consumer that parses coordinates out of
# gene names -- notably the official gapseqEval `predGeneEssentiality.GS`,
# which does sub("_[0-9]*:[0-9]*$","") and splits on ":" -- then sees a
# degenerate 1-bp range for every gene, so every knockout becomes a no-op and
# the model is scored as having zero essential genes.
#
# The rewrite itself is a deterministic, collision-free character
# substitution, hence invertible *if the original names are still available*.
# Every GEMAgent stage that runs purely in R (draft, gapfilled) saves its RDS
# with intact names, so such a file is the natural authority.
#
# Used by:
#   write_rds.R          -- repairing an RDS that was read back from SBML
#   export_sybil_model.R -- building a sybil modelorg for the evaluator
# ---------------------------------------------------------------------------

# The exact transformation the SBML writer applies to gene names.
sanitize_sbml_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)

# The *other* id mangling a round-trip can inflict. libSBML/cobra cannot put an
# illegal character in an `fbc:id`, so (unlike cobrar's writer, which substitutes
# "_") it spells the character out as "__<ord>__":
#
#     NZ_CP010052.1_410:1751  ->  NZ_CP010052__46__1_410__58__1751
#
# This form shows up when a model written by python-cobra is read back with
# sybilSBML::readSBMLmod, which does NOT undo the encoding. Which of the two
# substitutions is in play depends on which writer produced the SBML, so the
# lookup registers both keys and restoration works either way.
encode_sbml_id_cobra <- function(x) {
  x <- as.character(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    if (length(chars) == 0) return("")
    ok <- grepl("^[A-Za-z0-9_]$", chars)
    # !! utf8ToInt() 只接受**长度 1** 的字符串。把整个 chars 向量丢进去，
    #    它只取第一个元素并给出长度 1 的结果，sprintf() 再回收成整段长度，
    #    于是每个非法字符都被写成「第一个字符的码点」——例如
    #        NC_000964.3_3562565:3563542 -> NC_000964__78__3_3562565__78__3563542
    #    （78 是 "N"），而 python-cobra 实际写出的是 __46__ / __58__。
    #    key 因此永远命中不了，restore 静默退化为 no-op（坐标全丢、GE 全判非必需）。
    #    必须逐字符取码点。
    codes <- vapply(chars, utf8ToInt, integer(1), USE.NAMES = FALSE)
    paste(ifelse(ok, chars, sprintf("__%d__", codes)), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

# Inverse of encode_sbml_id_cobra(). Only plain printable ASCII is decoded, so a
# genuine name that happens to contain "__46__" is left untouched.
decode_sbml_id_cobra <- function(x) {
  x <- as.character(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    m <- gregexpr("__[0-9]+__", s)
    toks <- regmatches(s, m)[[1]]
    if (length(toks) == 0) return(s)
    regmatches(s, m) <- list(vapply(toks, function(t) {
      n <- suppressWarnings(as.integer(sub("^__([0-9]+)__$", "\\1", t)))
      if (is.na(n) || n < 33L || n > 126L) return(t)
      ch <- intToUtf8(n)
      if (grepl("^[A-Za-z0-9_]$", ch)) return(t)
      ch
    }, character(1), USE.NAMES = FALSE))
    s
  }, character(1), USE.NAMES = FALSE)
}

# ---- 第三种转义：cobrar::readSBMLmod 的 "__<ASCII 名称>__" ----------------
# libSBML 不允许 SId 里出现 "." / ":"，读入时要把它们改写掉。cobrar 的
# readSBMLmod 用的既不是 "_" 替换（cobrar 的 writer），也不是 "__<ord>__"
# （python-cobra 的 writer），而是**把字符名拼出来**并加 "gp_" 前缀：
#
#     NZ_CP010052.1_85737:86595
#       -> gp_NZ_CP__ZERO____ONE____ZERO____ZERO____FIVE____TWO____DOT____ONE_____EIGHT____FIVE____SEVEN____THREE____SEVEN____COLON____EIGHT____SIX____FIVE____NINE____FIVE__
#
# 这是**第三套**互不相通的写法，decode_sbml_id_cobra() 对它完全无效，于是
# seed_foreign_draft.R（用 cobrar::readSBMLmod 读外源 SBML）的解码会静默退化成
# no-op：draft 的基因名全是 "gp_...__DOT__..."，坐标解析率 0/855，
#   - GPR review 的 KOfam 桥（gpr_check_prompt.py --gene_pos_table）查不到基因，
#   - 评测端 postable 也对不上号 —— 整个外源 GEM 路线都会失效。
#
# COBRAR_CHAR_TOKEN 是实测枚举出来的：bsub 的 ModelSEED/CarveMe 坐标基因名里
# 只会出现下面 12 个 token（数字 + 点 + 冒号）。未列出的 "__XXX__" 一律原样保留，
# 以免把恰好长成这样的真名改坏。
COBRAR_CHAR_TOKEN <- c("0" = "ZERO", "1" = "ONE", "2" = "TWO", "3" = "THREE",
                       "4" = "FOUR", "5" = "FIVE", "6" = "SIX", "7" = "SEVEN",
                       "8" = "EIGHT", "9" = "NINE", "." = "DOT", ":" = "COLON")

encode_sbml_id_cobrar <- function(x) {
  x <- as.character(x)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    chars <- strsplit(s, "", fixed = TRUE)[[1]]
    if (length(chars) == 0) return("gp_")
    paste0("gp_", paste(ifelse(chars %in% names(COBRAR_CHAR_TOKEN),
                               paste0("__", COBRAR_CHAR_TOKEN[chars], "__"),
                               chars), collapse = ""))
  }, character(1), USE.NAMES = FALSE)
}

# 注意 token 的正则必须是 "__[A-Z]+__"（中间**不含下划线**）。用 "__[A-Z_]+__"
# 会因为贪婪匹配把 "__ZERO____ONE__" 当成一个大 token。
decode_sbml_id_cobrar <- function(x) {
  x <- as.character(x)
  rev_map <- stats::setNames(names(COBRAR_CHAR_TOKEN), COBRAR_CHAR_TOKEN)
  vapply(x, function(s) {
    if (is.na(s)) return(NA_character_)
    m <- gregexpr("__[A-Z]+__", s)
    toks <- regmatches(s, m)[[1]]
    if (length(toks)) {
      regmatches(s, m) <- list(vapply(toks, function(t) {
        key <- sub("^__(.*)__$", "\\1", t)
        if (key %in% names(rev_map)) rev_map[[key]] else t
      }, character(1), USE.NAMES = FALSE))
    }
    sub("^gp_", "", s)
  }, character(1), USE.NAMES = FALSE)
}

# 一站式解码：一个名字只会是三种转义中的一种，按顺序试即可（三者互不干扰）。
decode_sbml_id_any <- function(x) decode_sbml_id_cobrar(decode_sbml_id_cobra(x))

# 自检：encode 与 decode 必须严格互逆，否则 gene_name_lookup 的 "__<ord>__"
# key 一定对不上 python-cobra 写出的 SBML，restore 会静默退化成 no-op。
# 放在这里（两个函数都定义完之后）而不是文件顶端，避免 source 时前向引用。
.stop_if_encode_is_broken <- function(probe = c("NZ_CP010052.1_410:1751",
                                                "NC_000964.3_3562565:3563542")) {
  enc <- encode_sbml_id_cobra(probe)
  dec <- decode_sbml_id_cobra(enc)
  if (!identical(dec, probe)) {
    stop(sprintf(paste0("restore_gene_names.R: encode/decode 不是互逆的，",
                        "SBML id 转义将无法还原。\n  in : %s\n  enc: %s\n  out: %s"),
                 paste(probe, collapse = ", "), paste(enc, collapse = ", "),
                 paste(dec, collapse = ", ")))
  }
  if (!identical(enc[1], "NZ_CP010052__46__1_410__58__1751")) {
    stop("restore_gene_names.R: encode_sbml_id_cobra 的输出与 python-cobra 的 ",
         "SBML id 转义不一致（期望 __46__ / __58__）。")
  }
  # 第三套转义（cobrar::readSBMLmod）同样必须互逆，否则用 cobrar 读外源 SBML 时
  #（seed_foreign_draft.R）坐标会全丢，GPR review 与 GE 评测一并失效。
  enc.cobrar <- encode_sbml_id_cobrar(probe)
  dec.cobrar <- decode_sbml_id_cobrar(enc.cobrar)
  if (!identical(dec.cobrar, probe)) {
    stop(sprintf(paste0("restore_gene_names.R: cobrar 转义不是互逆的。\n",
                        "  in : %s\n  enc: %s\n  out: %s"),
                 paste(probe, collapse = ", "), paste(enc.cobrar, collapse = ", "),
                 paste(dec.cobrar, collapse = ", ")))
  }
  invisible(TRUE)
}
.stop_if_encode_is_broken()

# Build `mangled name -> true name` from any vector of intact names.
# Guards against two distinct true names collapsing onto one key, which would
# make the inverse ambiguous; the first occurrence wins, the rest are dropped.
#
# Six keys are registered for every name: the bare and "G_"-prefixed forms of
# the "_" sanitized id, the "__<ord>__" cobra id, and the "gp_...__NAME__"
# cobrar id. cobrar::readSBMLmod strips the "G_" prefix while sybilSBML::readSBMLmod
# keeps it, and any of the three writer substitutions may appear, so registering
# all six lets the same lookup serve every reader/writer combination.
gene_name_lookup <- function(auth.names) {
  auth.names <- as.character(auth.names)
  keys <- c(sanitize_sbml_id(auth.names), encode_sbml_id_cobra(auth.names),
            encode_sbml_id_cobrar(auth.names))
  vals <- c(auth.names, auth.names, auth.names)
  keys <- c(keys, paste0("G_", keys))
  vals <- c(vals, vals)
  if (anyDuplicated(keys)) {
    n.dup <- sum(duplicated(keys))
    warning(sprintf(paste0("%d gene-name key(s) collide after SBML id mangle; ",
                           "keeping the first occurrence and leaving the rest ",
                           "unmapped."), n.dup))
    keep <- !duplicated(keys)
    keys <- keys[keep]
    vals <- vals[keep]
  }
  stats::setNames(vals, keys)
}
# Convenience: same thing, read from an RDS whose gene names are intact.
# Falls back through the candidate list so callers can hand over several.
gene_name_lookup_from_rds <- function(paths) {
  for (p in paths) {
    if (is.character(p) && length(p) == 1 && nzchar(p) && file.exists(p)) {
      lookup <- gene_name_lookup(readRDS(p)@allGenes)
      attr(lookup, "source") <- p
      return(lookup)
    }
  }
  NULL
}

# Rewrite @allGenes and @genes IN PLACE (order-preserving, so the index-based
# @gprRules stay valid). Both slots must be rewritten together: sybil resolves
# `optimizeProb(gene = ...)` through @allGenes, but builds the knockout from
# @genes, so touching only one of them silently disables the knockout.
restore_gene_names <- function(mod, lookup, verbose = FALSE) {
  remap <- function(v) {
    if (length(v) == 0) return(v)
    hit <- lookup[v]
    ifelse(is.na(hit), v, unname(hit))
  }
  n.before <- length(mod@allGenes)
  mod@allGenes <- as.character(remap(mod@allGenes))
  mod@genes    <- lapply(mod@genes, remap)

  if (verbose) {
    still <- is.na(lookup[sanitize_sbml_id(mod@allGenes)])
    cat(sprintf("  gene names: %d/%d restored%s\n",
                n.before - sum(still), n.before,
                if (!is.null(attr(lookup, "source")))
                  paste0(" from ", basename(attr(lookup, "source"))) else ""))
    if (any(still)) {
      cat(sprintf("    still sanitized (%d): %s%s\n", sum(still),
                  paste(head(mod@allGenes[still], 5), collapse = ", "),
                  if (sum(still) > 5) " ..." else ""))
    }
    rg <- unique(unlist(mod@genes))
    rg <- rg[!is.na(rg) & nzchar(rg)]
    bad <- setdiff(rg, mod@allGenes)
    cat(sprintf("    self-check: @genes entries missing from @allGenes = %d%s\n",
                length(bad),
                if (length(bad)) paste0(" (e.g. ", paste(head(bad, 3), collapse = ", "), ")") else ""))
  }
  mod
}