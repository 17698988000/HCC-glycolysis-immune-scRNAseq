# ============================================================
# scripts/utils/locked_fig2c_groups.R
#
# Authoritative Figure 2C group assignment helper.
#
# The submission package uses a rank-balanced split of the final
# 15,391 tumor-derived hepatocytes. Two cells share the median AUCell
# score, so downstream scripts MUST read the locked assignment table
# instead of recreating groups with a simple > median comparison.
# ============================================================

LOCKED_FIG2C_MEDIAN_AUC <- 0.212617779598525
LOCKED_FIG2C_TUMOR_HEPATOCYTE_N <- 15391L
LOCKED_FIG2C_GLYCOHIGH_N <- 7695L
LOCKED_FIG2C_GLYCOLOW_N <- 7696L
LOCKED_FIG2C_MEDIAN_TIE_N <- 2L

find_repository_root <- function(start = getwd()) {
  start <- normalizePath(start, winslash = "/", mustWork = TRUE)
  candidates <- unique(c(
    start,
    dirname(start),
    dirname(dirname(start)),
    Sys.getenv("PROJECT_DIR", unset = "")
  ))
  candidates <- candidates[nzchar(candidates)]
  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "README.md")) &&
        dir.exists(file.path(candidate, "locked_source_data"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Could not locate repository root. Run from the repository root or set PROJECT_DIR.",
    call. = FALSE
  )
}

load_locked_fig2c_assignments <- function(repository_root = find_repository_root()) {
  path <- file.path(
    repository_root,
    "locked_source_data",
    "single_cell",
    "Fig2C_15391_tumor_hepatocyte_GlycoHigh_GlycoLow_FIXED.csv"
  )
  if (!file.exists(path)) stop("Locked Figure 2C assignment file not found: ", path, call. = FALSE)
  x <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("cell", "Glycolysis_AUC", "GlycoGroup")
  missing <- setdiff(required, names(x))
  if (length(missing) > 0) stop("Locked Figure 2C file is missing: ", paste(missing, collapse = ", "), call. = FALSE)
  x$cell <- as.character(x$cell)
  x$Glycolysis_AUC <- as.numeric(x$Glycolysis_AUC)
  x$GlycoGroup <- as.character(x$GlycoGroup)
  if (nrow(x) != LOCKED_FIG2C_TUMOR_HEPATOCYTE_N) stop("Locked Figure 2C row count mismatch.", call. = FALSE)
  counts <- table(x$GlycoGroup)
  if (!identical(as.integer(counts[["GlycoHigh"]]), LOCKED_FIG2C_GLYCOHIGH_N) ||
      !identical(as.integer(counts[["GlycoLow"]]), LOCKED_FIG2C_GLYCOLOW_N)) {
    stop("Locked Figure 2C group-count mismatch.", call. = FALSE)
  }
  observed_median <- stats::median(x$Glycolysis_AUC, na.rm = TRUE)
  if (abs(observed_median - LOCKED_FIG2C_MEDIAN_AUC) > 1e-12) stop("Locked Figure 2C median mismatch.", call. = FALSE)
  tie_n <- sum(abs(x$Glycolysis_AUC - LOCKED_FIG2C_MEDIAN_AUC) < 1e-12, na.rm = TRUE)
  if (tie_n != LOCKED_FIG2C_MEDIAN_TIE_N) stop("Locked Figure 2C median-tie count mismatch.", call. = FALSE)
  x
}

map_locked_fig2c_groups <- function(cell_ids, high_label = "High", low_label = "Low", repository_root = find_repository_root()) {
  locked <- load_locked_fig2c_assignments(repository_root)
  group <- locked$GlycoGroup[match(as.character(cell_ids), locked$cell)]
  out <- ifelse(group == "GlycoHigh", high_label, ifelse(group == "GlycoLow", low_label, NA_character_))
  out
}
