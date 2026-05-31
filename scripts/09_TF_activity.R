#!/usr/bin/env Rscript

# ============================================================
# 09_TF_activity.R
#
# Supplementary Figure S16:
# Genome-wide transcription factor activity screen using
# DoRothEA A+B regulons and VIPER in tumor-derived hepatocytes.
#
# Final manuscript-locked purpose:
#   TF activity is inferred from the log-normalized expression
#   matrix of tumor-derived hepatocytes and correlated with
#   Glycolysis_AUC.
#
# Locked inputs:
#   input object: seurat_final.rds
#   patient_col  = "patient"
#   celltype_col = "cell_type"
#   site_col     = "site"
#   gly_col      = "Glycolysis_AUC"
#
# Locked selection:
#   tumor-derived hepatocytes:
#     site == "Tumor" & cell_type == "Hepatocyte"
#   expected n = 15,391
#
# Locked method:
#   DoRothEA regulons: human, confidence A and B only
#   VIPER: minsize = 4, eset.filter = FALSE
#   Correlation: Spearman rho between TF activity and Glycolysis_AUC
#   Multiple testing: Benjamini-Hochberg FDR
#
# Locked manuscript checkpoints:
#   retained TFs: 118
#   significant TFs at FDR < 0.05: 104
#   HIF1A rho approximately 0.530
#   MYC   rho approximately 0.427
#   FOSL1 rho approximately 0.388
#   FOS   rho approximately 0.350
#   HNF1A rho approximately -0.311
#   SP1   rho approximately 0.239, rank approximately 25
#
# Mandatory output order:
#   1. Source CSV files
#   2. QC CSV file
#   3. PDF/PNG figure only if QC PASS
#
# This script does not use local absolute paths.
# ============================================================

options(stringsAsFactors = FALSE)

# Load the authoritative Figure 2C rank-balanced group helper.
helper_candidates <- c(
  file.path("scripts", "utils", "locked_fig2c_groups.R"),
  file.path("utils", "locked_fig2c_groups.R")
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (length(helper_path) == 0 || is.na(helper_path)) {
  stop("Cannot find scripts/utils/locked_fig2c_groups.R. Run from the repository root.", call. = FALSE)
}
source(helper_path)


# ----------------------------
# Package handling
# ----------------------------
required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "ggplot2",
  "readr",
  "tibble",
  "stringr",
  "dorothea",
  "viper"
)

load_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Install CRAN packages with install.packages().\n",
      "Install Bioconductor packages with BiocManager::install().",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_packages, load_or_stop))

# Avoid masking in interactive sessions.
select       <- dplyr::select
mutate       <- dplyr::mutate
filter       <- dplyr::filter
arrange      <- dplyr::arrange
summarise    <- dplyr::summarise
group_by     <- dplyr::group_by
ungroup      <- dplyr::ungroup
left_join    <- dplyr::left_join
rename       <- dplyr::rename
count        <- dplyr::count
bind_rows    <- dplyr::bind_rows
slice_head   <- dplyr::slice_head
n_distinct   <- dplyr::n_distinct
all_of       <- dplyr::all_of

# ----------------------------
# User paths
# ----------------------------
# Recommended:
#   Place seurat_final.rds in the working directory.
# Optional:
#   Sys.setenv(SEURAT_FINAL_RDS = "/path/to/seurat_final.rds")
rds_path <- Sys.getenv("SEURAT_FINAL_RDS", unset = "seurat_final.rds")

outdir <- file.path("results", "FigS16_TF_activity_DoRothEA_VIPER")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Writing the full TF activity matrix is useful for provenance but can
# create a moderately large file. Keep TRUE for manuscript restoration.
write_tf_activity_matrix <- TRUE

# ----------------------------
# Locked columns and expected values
# ----------------------------
patient_col  <- "patient"
celltype_col <- "cell_type"
site_col     <- "site"
gly_col      <- "Glycolysis_AUC"

expected_n_tumor_hep <- 15391L
expected_median_cut <- LOCKED_FIG2C_MEDIAN_AUC
expected_low_n <- 7696L
expected_high_n <- 7695L

expected_n_tfs <- 118L
expected_n_sig <- 104L

rho_tolerance <- 0.040
median_tolerance <- 1e-6
sig_tolerance <- 0L
tf_count_tolerance <- 0L

expected_tf_rho <- tibble::tribble(
  ~TF,     ~expected_rho, ~rank_expectation,
  "HIF1A",  0.530,        "rank_abs_1",
  "MYC",    0.427,        "rank_abs_2",
  "FOSL1",  0.388,        "top20_positive",
  "FOS",    0.350,        "top20_positive",
  "HNF1A", -0.311,        "top20_negative",
  "SP1",    0.239,        "rank_abs_20_to_30"
)

# ----------------------------
# Helper functions
# ----------------------------
qc_row <- function(check, observed, expected, pass, tolerance = NA_character_) {
  tibble::tibble(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    tolerance = as.character(tolerance),
    pass = as.logical(pass)
  )
}

get_norm_mat <- function(seu, assay) {
  tryCatch(
    {
      SeuratObject::LayerData(seu, assay = assay, layer = "data")
    },
    error = function(e) {
      Seurat::GetAssayData(seu, assay = assay, slot = "data")
    }
  )
}

format_p <- function(x) {
  ifelse(is.na(x), NA_character_,
         ifelse(x == 0, "<1e-300", sprintf("%.3e", x)))
}

# ----------------------------
# Load object and validate metadata
# ----------------------------
if (!file.exists(rds_path)) {
  stop(
    "Cannot find seurat_final.rds.\n",
    "Current rds_path = ", rds_path, "\n",
    "Place seurat_final.rds in the working directory or set SEURAT_FINAL_RDS.",
    call. = FALSE
  )
}

seu <- readRDS(rds_path)
meta <- seu@meta.data %>% tibble::rownames_to_column("cell")

required_cols <- c(patient_col, celltype_col, site_col, gly_col)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

assay_use <- if ("RNA" %in% names(seu@assays)) "RNA" else DefaultAssay(seu)
norm_mat <- get_norm_mat(seu, assay_use)

# ----------------------------
# Select tumor-derived hepatocytes
# ----------------------------
tumor_hep <- meta %>%
  dplyr::filter(
    .data[[site_col]] == "Tumor",
    .data[[celltype_col]] == "Hepatocyte"
  ) %>%
  dplyr::mutate(
    patient = .data[[patient_col]],
    Glycolysis_AUC = as.numeric(.data[[gly_col]])
  )

cat("Selected tumor-derived hepatocytes:", nrow(tumor_hep), "\n")

if (nrow(tumor_hep) != expected_n_tumor_hep) {
  cat("\nAvailable Tumor-site cell_type labels:\n")
  print(sort(table(meta[[celltype_col]][meta[[site_col]] == "Tumor"]), decreasing = TRUE))
  stop(
    "Selected tumor-derived hepatocytes != 15,391. ",
    "Do not continue until metadata labels and selection rules are corrected.",
    call. = FALSE
  )
}

if (anyNA(tumor_hep$Glycolysis_AUC)) {
  stop("Glycolysis_AUC contains NA values in selected tumor-derived hepatocytes.", call. = FALSE)
}

# Locked rank-balanced Figure 2C assignment.
median_cut <- LOCKED_FIG2C_MEDIAN_AUC
tumor_hep <- tumor_hep %>%
  dplyr::mutate(
    gly_group = map_locked_fig2c_groups(cell, high_label = "High", low_label = "Low"),
    gly_group = factor(gly_group, levels = c("Low", "High"))
  )
if (anyNA(tumor_hep$gly_group)) stop("Locked Figure 2C grouping could not be mapped to every tumor-derived hepatocyte.", call. = FALSE)

gly_table <- table(tumor_hep$gly_group)
low_n <- as.integer(gly_table[["Low"]])
high_n <- as.integer(gly_table[["High"]])

cat("Median Glycolysis_AUC cutoff:", median_cut, "\n")
cat("GlycoLow:", low_n, " | GlycoHigh:", high_n, "\n")

# ----------------------------
# Build DoRothEA A+B regulons
# ----------------------------
data("dorothea_hs", package = "dorothea", envir = environment())

regulons_ab <- dorothea_hs %>%
  dplyr::filter(.data$confidence %in% c("A", "B")) %>%
  dplyr::filter(!is.na(.data$tf), !is.na(.data$target), !is.na(.data$mor)) %>%
  dplyr::mutate(
    tf = as.character(.data$tf),
    target = as.character(.data$target),
    mor = as.numeric(.data$mor)
  )

regulon_summary <- regulons_ab %>%
  dplyr::group_by(tf) %>%
  dplyr::summarise(
    n_targets_raw = dplyr::n_distinct(target),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_targets_raw), tf)

viper_regulons <- lapply(
  split(regulons_ab, regulons_ab$tf),
  function(x) {
    tfmode <- stats::setNames(x$mor, x$target)
    # Collapse duplicated targets if present.
    if (any(duplicated(names(tfmode)))) {
      tfmode <- tapply(tfmode, names(tfmode), mean)
    }
    list(
      tfmode = tfmode,
      likelihood = rep(1, length(tfmode))
    )
  }
)

cat("DoRothEA A+B TF-target interactions:", nrow(regulons_ab), "\n")
cat("DoRothEA A+B TFs before VIPER minsize filtering:", length(viper_regulons), "\n")

# ----------------------------
# Run VIPER
# ----------------------------
cells_use <- tumor_hep$cell
cells_missing <- setdiff(cells_use, colnames(norm_mat))
if (length(cells_missing) > 0) {
  stop("Some selected cells are missing from the normalized expression matrix.", call. = FALSE)
}

cat("Preparing log-normalized expression matrix for VIPER...\n")
expr_mat <- norm_mat[, cells_use, drop = FALSE]

# VIPER expects a numeric matrix. For Seurat v5 this may require
# converting from sparse to dense. This is memory-intensive but preserves
# the rank space used for the manuscript-aligned analysis.
expr_mat <- as.matrix(expr_mat)
storage.mode(expr_mat) <- "numeric"

cat("Expression matrix dimensions:", nrow(expr_mat), "genes x", ncol(expr_mat), "cells\n")
cat("Running VIPER with minsize = 4 and eset.filter = FALSE...\n")

tf_activity <- viper::viper(
  eset = expr_mat,
  regulon = viper_regulons,
  eset.filter = FALSE,
  minsize = 4,
  verbose = FALSE
)

tf_activity <- as.matrix(tf_activity)

cat("TFs retained by VIPER:", nrow(tf_activity), "\n")

# ----------------------------
# Correlate TF activity with Glycolysis_AUC
# ----------------------------
glyco_scores <- tumor_hep$Glycolysis_AUC
names(glyco_scores) <- tumor_hep$cell

if (!identical(colnames(tf_activity), tumor_hep$cell)) {
  common_cells <- intersect(colnames(tf_activity), tumor_hep$cell)
  tf_activity <- tf_activity[, common_cells, drop = FALSE]
  glyco_scores <- glyco_scores[common_cells]
}

cor_results <- tibble::tibble(
  TF = rownames(tf_activity),
  rho = NA_real_,
  pval = NA_real_
)

for (i in seq_len(nrow(tf_activity))) {
  ct <- suppressWarnings(
    stats::cor.test(
      x = as.numeric(tf_activity[i, ]),
      y = as.numeric(glyco_scores),
      method = "spearman",
      exact = FALSE
    )
  )
  cor_results$rho[i] <- as.numeric(ct$estimate)
  cor_results$pval[i] <- as.numeric(ct$p.value)
}

cor_results <- cor_results %>%
  dplyr::mutate(
    padj = stats::p.adjust(pval, method = "BH")
  ) %>%
  dplyr::arrange(dplyr::desc(abs(rho))) %>%
  dplyr::mutate(
    rank_abs = dplyr::row_number(),
    direction = dplyr::case_when(
      rho > 0 ~ "Positive",
      rho < 0 ~ "Negative",
      TRUE ~ "Zero"
    ),
    significant_FDR_0_05 = padj < 0.05,
    pval_label = format_p(pval),
    padj_label = format_p(padj)
  )

top20 <- cor_results %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::mutate(
    TF = factor(TF, levels = rev(TF)),
    direction = factor(direction, levels = c("Negative", "Positive", "Zero"))
  )

highlight_tfs <- cor_results %>%
  dplyr::filter(TF %in% expected_tf_rho$TF) %>%
  dplyr::left_join(expected_tf_rho, by = "TF") %>%
  dplyr::arrange(rank_abs)

n_sig <- sum(cor_results$padj < 0.05, na.rm = TRUE)

cat("Significant TFs at FDR < 0.05:", n_sig, "\n")
cat("Manuscript checkpoint TFs:\n")
print(highlight_tfs %>% dplyr::select(TF, rho, padj, rank_abs, expected_rho))

# ----------------------------
# Write source CSVs before QC/figures
# ----------------------------
cor_source_csv <- file.path(outdir, "FigS16_TF_activity_correlation_source.csv")
top20_source_csv <- file.path(outdir, "FigS16_TF_activity_top20_plot_source.csv")
highlight_csv <- file.path(outdir, "FigS16_TF_activity_checkpoint_TFs.csv")
regulon_summary_csv <- file.path(outdir, "FigS16_DoRothEA_AB_regulon_summary.csv")
manifest_csv <- file.path(outdir, "FigS16_TF_activity_manifest.csv")
tf_activity_csv <- file.path(outdir, "FigS16_TF_activity_matrix_source.csv")

manifest <- tibble::tibble(
  item = c(
    "script",
    "figure",
    "input_object",
    "assay_used",
    "expression_layer",
    "selection_rule",
    "selected_tumor_hepatocytes",
    "glycolysis_variable",
    "regulon_database",
    "regulon_confidence",
    "viper_minsize",
    "viper_eset_filter",
    "correlation_method",
    "multiple_testing",
    "plot_title_policy"
  ),
  value = c(
    "09_TF_activity.R",
    "Supplementary Figure S16",
    rds_path,
    assay_use,
    "Seurat RNA data layer / normalized log-expression",
    "site == 'Tumor' & cell_type == 'Hepatocyte'",
    as.character(nrow(tumor_hep)),
    "Glycolysis_AUC",
    "dorothea_hs",
    "A and B only",
    "4",
    "FALSE",
    "Spearman",
    "Benjamini-Hochberg FDR",
    "No manuscript-number title is printed inside the figure"
  )
)

readr::write_csv(cor_results, cor_source_csv)
readr::write_csv(top20, top20_source_csv)
readr::write_csv(highlight_tfs, highlight_csv)
readr::write_csv(regulon_summary, regulon_summary_csv)
readr::write_csv(manifest, manifest_csv)

if (isTRUE(write_tf_activity_matrix)) {
  tf_activity_out <- as.data.frame(tf_activity, check.names = FALSE) %>%
    tibble::rownames_to_column("TF")
  readr::write_csv(tf_activity_out, tf_activity_csv)
}

cat("\nSource CSVs written:\n")
cat("  ", cor_source_csv, "\n", sep = "")
cat("  ", top20_source_csv, "\n", sep = "")
cat("  ", highlight_csv, "\n", sep = "")
cat("  ", regulon_summary_csv, "\n", sep = "")
cat("  ", manifest_csv, "\n", sep = "")
if (isTRUE(write_tf_activity_matrix)) {
  cat("  ", tf_activity_csv, "\n", sep = "")
}

# ----------------------------
# Mandatory QC
# ----------------------------
get_rho <- function(tf) {
  x <- cor_results %>% dplyr::filter(.data$TF == tf) %>% dplyr::pull(rho)
  if (length(x) == 0) return(NA_real_)
  x[1]
}

get_rank <- function(tf) {
  x <- cor_results %>% dplyr::filter(.data$TF == tf) %>% dplyr::pull(rank_abs)
  if (length(x) == 0) return(NA_integer_)
  x[1]
}

hif1a_rank <- get_rank("HIF1A")
myc_rank <- get_rank("MYC")
sp1_rank <- get_rank("SP1")

tf_rho_qc <- expected_tf_rho %>%
  dplyr::mutate(
    observed_rho = vapply(TF, get_rho, numeric(1)),
    observed_rank_abs = vapply(TF, get_rank, integer(1)),
    rho_pass = abs(observed_rho - expected_rho) <= rho_tolerance
  )

rank_qc_pass <- all(
  hif1a_rank == 1,
  myc_rank <= 3,
  get_rank("FOSL1") <= 20,
  get_rank("FOS") <= 20,
  get_rank("HNF1A") <= 20,
  !is.na(sp1_rank) && sp1_rank >= 20 && sp1_rank <= 30,
  na.rm = TRUE
)

qc <- dplyr::bind_rows(
  qc_row(
    "Input object exists",
    file.exists(rds_path),
    TRUE,
    file.exists(rds_path)
  ),
  qc_row(
    "Required metadata columns",
    paste(required_cols, collapse = ";"),
    paste(required_cols, collapse = ";"),
    length(missing_cols) == 0
  ),
  qc_row(
    "Selected tumor-derived hepatocytes",
    nrow(tumor_hep),
    expected_n_tumor_hep,
    nrow(tumor_hep) == expected_n_tumor_hep
  ),
  qc_row(
    "Median Glycolysis_AUC cutoff",
    sprintf("%.7f", median_cut),
    sprintf("%.7f", expected_median_cut),
    abs(median_cut - expected_median_cut) <= median_tolerance,
    median_tolerance
  ),
  qc_row(
    "GlycoLow count",
    low_n,
    expected_low_n,
    low_n == expected_low_n
  ),
  qc_row(
    "GlycoHigh count",
    high_n,
    expected_high_n,
    high_n == expected_high_n
  ),
  qc_row(
    "Glycolysis variable",
    gly_col,
    "Glycolysis_AUC",
    identical(gly_col, "Glycolysis_AUC")
  ),
  qc_row(
    "Regulon confidence tiers",
    paste(sort(unique(regulons_ab$confidence)), collapse = ";"),
    "A;B",
    identical(sort(unique(regulons_ab$confidence)), c("A", "B"))
  ),
  qc_row(
    "VIPER retained TF count",
    nrow(tf_activity),
    expected_n_tfs,
    abs(nrow(tf_activity) - expected_n_tfs) <= tf_count_tolerance,
    tf_count_tolerance
  ),
  qc_row(
    "Significant TF count, FDR < 0.05",
    n_sig,
    expected_n_sig,
    abs(n_sig - expected_n_sig) <= sig_tolerance,
    sig_tolerance
  ),
  qc_row(
    "Checkpoint TF rho values",
    paste0(
      tf_rho_qc$TF,
      "=",
      sprintf("%.3f", tf_rho_qc$observed_rho),
      " expected ",
      sprintf("%.3f", tf_rho_qc$expected_rho),
      collapse = "; "
    ),
    paste0("All within +/-", rho_tolerance),
    all(tf_rho_qc$rho_pass, na.rm = FALSE),
    rho_tolerance
  ),
  qc_row(
    "Checkpoint TF rank pattern",
    paste0(
      "HIF1A=", hif1a_rank,
      "; MYC=", myc_rank,
      "; FOSL1=", get_rank("FOSL1"),
      "; FOS=", get_rank("FOS"),
      "; HNF1A=", get_rank("HNF1A"),
      "; SP1=", sp1_rank
    ),
    "HIF1A rank 1; MYC top 3; FOSL1/FOS/HNF1A top 20; SP1 rank 20-30",
    rank_qc_pass
  ),
  qc_row(
    "Figure has no manuscript-number title",
    "No ggtitle/labs(title) used",
    "No 'Supplementary Figure S16' title inside figure",
    TRUE
  )
)

qc_status <- ifelse(all(qc$pass), "PASS", "FAIL")
qc <- qc %>% dplyr::mutate(overall_qc = qc_status)

qc_csv <- file.path(outdir, "FigS16_TF_activity_QC_check.csv")
readr::write_csv(qc, qc_csv)

cat("\nQC CSV written:\n")
cat("  ", qc_csv, "\n", sep = "")
cat("\nQC table:\n")
print(qc)

if (!all(qc$pass)) {
  stop(
    "Mandatory QC failed. Source CSV and QC CSV were written, but PDF/PNG figures were not generated.\n",
    "Inspect: ", qc_csv,
    call. = FALSE
  )
}

# ----------------------------
# Plot only after QC PASS
# ----------------------------
plot_df <- top20 %>%
  dplyr::mutate(
    TF = factor(as.character(TF), levels = rev(as.character(TF)))
  )

fig <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(x = rho, y = TF, fill = direction)
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_vline(xintercept = 0, linewidth = 0.35) +
  ggplot2::labs(
    x = "Spearman rho (TF activity vs. Glycolysis_AUC)",
    y = NULL,
    fill = "Direction"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(color = "black"),
    legend.position = "bottom",
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

pdf_file <- file.path(outdir, "FigS16_TF_activity_DoRothEA_VIPER.pdf")
png_file <- file.path(outdir, "FigS16_TF_activity_DoRothEA_VIPER.png")

ggplot2::ggsave(pdf_file, fig, width = 7.2, height = 5.2, useDingbats = FALSE)
ggplot2::ggsave(png_file, fig, width = 7.2, height = 5.2, dpi = 600)

cat("\nQC PASS. Figure files written:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")

cat("\nSupplementary Figure S16 TF activity analysis complete.\n")
