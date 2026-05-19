#!/usr/bin/env Rscript

# ============================================================
# 08_glycolysis_gradient.R
#
# Supplementary Figure S15 glycolysis-gradient analysis.
#
# Final manuscript-locked purpose:
#   Tumor-derived hepatocytes (site == "Tumor" &
#   cell_type == "Hepatocyte"; n = 15,391) are stratified into
#   10 equal-width bins spanning Glycolysis_AUC. For each bin,
#   the mean log-normalized expression of ENO1, LDHA, SPP1, MIF,
#   and PTGES is computed and plotted.
#
# Locked inputs:
#   input object: seurat_final.rds
#   patient_col  = "patient"
#   celltype_col = "cell_type"
#   site_col     = "site"
#   gly_col      = "Glycolysis_AUC"
#
# Locked tumor-hepatocyte stratification checks:
#   selected tumor-derived hepatocytes = 15,391
#   median Glycolysis_AUC cutoff       = 0.2203849
#   GlycoLow  = 7,696
#   GlycoHigh = 7,695
#   High rule = Glycolysis_AUC > median_cut
#
# Mandatory output order:
#   1. Source CSVs
#   2. QC CSV
#   3. PDF/PNG only if QC PASS
#
# This script intentionally does NOT use ntile() quantile bins.
# S15 uses 10 equal-width AUCell-score bins.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
  library(scales)
})

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
slice        <- dplyr::slice
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

outdir <- file.path("results", "FigS15_glycolysis_gradient")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Locked columns and targets
# ----------------------------
patient_col  <- "patient"
celltype_col <- "cell_type"
site_col     <- "site"
gly_col      <- "Glycolysis_AUC"

target_genes <- c("ENO1", "LDHA", "SPP1", "MIF", "PTGES")
glycolytic_genes <- c("ENO1", "LDHA")
ligand_genes <- c("SPP1", "MIF", "PTGES")

expected_n_tumor_hep <- 15391L
expected_low_n <- 7696L
expected_high_n <- 7695L
expected_median_cut <- 0.2203849
median_tolerance <- 1e-6

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

pull_gene_expr <- function(genes, cells, norm_mat) {
  genes_present <- intersect(genes, rownames(norm_mat))

  if (length(genes_present) != length(genes)) {
    warning("Missing genes: ", paste(setdiff(genes, genes_present), collapse = ", "))
  }

  if (length(genes_present) == 0) {
    stop("None of the target genes were found in the normalized expression matrix.", call. = FALSE)
  }

  expr <- as.matrix(norm_mat[genes_present, cells, drop = FALSE])
  expr <- t(expr)
  expr <- as.data.frame(expr)
  expr <- tibble::rownames_to_column(expr, "cell")
  expr
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

missing_target_genes <- setdiff(target_genes, rownames(norm_mat))
if (length(missing_target_genes) > 0) {
  stop("Missing target genes in normalized expression matrix: ",
       paste(missing_target_genes, collapse = ", "), call. = FALSE)
}

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
    "Do not continue plotting until metadata labels and selection rules are corrected.",
    call. = FALSE
  )
}

if (anyNA(tumor_hep$Glycolysis_AUC)) {
  stop("Glycolysis_AUC contains NA values in selected tumor-derived hepatocytes.", call. = FALSE)
}

# ----------------------------
# Locked median split check
# ----------------------------
median_cut <- median(tumor_hep$Glycolysis_AUC, na.rm = TRUE)

tumor_hep <- tumor_hep %>%
  dplyr::mutate(
    gly_group = ifelse(Glycolysis_AUC > median_cut, "High", "Low"),
    gly_group = factor(gly_group, levels = c("Low", "High"))
  )

gly_table <- table(tumor_hep$gly_group)
low_n <- as.integer(gly_table[["Low"]])
high_n <- as.integer(gly_table[["High"]])

cat("Median Glycolysis_AUC cutoff:", median_cut, "\n")
cat("GlycoLow:", low_n, " | GlycoHigh:", high_n, "\n")

# ----------------------------
# Equal-width AUCell-score bins
# ----------------------------
gly_min <- min(tumor_hep$Glycolysis_AUC, na.rm = TRUE)
gly_max <- max(tumor_hep$Glycolysis_AUC, na.rm = TRUE)

if (!is.finite(gly_min) || !is.finite(gly_max) || gly_min == gly_max) {
  stop("Invalid Glycolysis_AUC range for equal-width binning.", call. = FALSE)
}

bin_breaks <- seq(gly_min, gly_max, length.out = 11)

tumor_hep <- tumor_hep %>%
  dplyr::mutate(
    glyco_bin = cut(
      Glycolysis_AUC,
      breaks = bin_breaks,
      include.lowest = TRUE,
      right = TRUE,
      labels = paste0("Bin", 1:10)
    ),
    glyco_bin_index = as.integer(glyco_bin),
    bin_label = paste0("Bin ", glyco_bin_index)
  )

bin_counts <- tumor_hep %>%
  dplyr::count(glyco_bin_index, bin_label, name = "n_cells") %>%
  dplyr::arrange(glyco_bin_index)

cat("\nEqual-width bin counts:\n")
print(bin_counts)

# ----------------------------
# Pull expression and write source data
# ----------------------------
expr_df <- pull_gene_expr(target_genes, tumor_hep$cell, norm_mat)

source_wide <- tumor_hep %>%
  dplyr::select(
    cell,
    patient,
    Glycolysis_AUC,
    gly_group,
    glyco_bin_index,
    bin_label
  ) %>%
  dplyr::left_join(expr_df, by = "cell")

source_long <- source_wide %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(target_genes),
    names_to = "gene",
    values_to = "log_normalized_expression"
  ) %>%
  dplyr::mutate(
    gene_class = dplyr::case_when(
      gene %in% glycolytic_genes ~ "Glycolytic enzyme",
      gene %in% ligand_genes ~ "Immunosuppressive ligand",
      TRUE ~ "Other"
    )
  )

bin_means <- source_long %>%
  dplyr::group_by(glyco_bin_index, bin_label, gene, gene_class) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_log_normalized_expression = mean(log_normalized_expression, na.rm = TRUE),
    median_log_normalized_expression = median(log_normalized_expression, na.rm = TRUE),
    sd_log_normalized_expression = stats::sd(log_normalized_expression, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(gene_class, gene, glyco_bin_index)

bin_ranges <- tumor_hep %>%
  dplyr::group_by(glyco_bin_index, bin_label) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    auc_min = min(Glycolysis_AUC, na.rm = TRUE),
    auc_max = max(Glycolysis_AUC, na.rm = TRUE),
    auc_mean = mean(Glycolysis_AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(glyco_bin_index)

manifest <- tibble::tibble(
  item = c(
    "script",
    "figure",
    "input_object",
    "assay_used",
    "selection_rule",
    "selected_tumor_hepatocytes",
    "binning_rule",
    "genes",
    "expression_layer",
    "plot_title_policy"
  ),
  value = c(
    "08_glycolysis_gradient.R",
    "Supplementary Figure S15",
    rds_path,
    assay_use,
    "site == 'Tumor' & cell_type == 'Hepatocyte'",
    as.character(nrow(tumor_hep)),
    "10 equal-width bins spanning Glycolysis_AUC range in selected tumor-derived hepatocytes",
    paste(target_genes, collapse = ";"),
    "Seurat RNA data layer / normalized expression",
    "No 'Supplementary Figure S15' title is printed inside the figure"
  )
)

cell_source_csv <- file.path(outdir, "FigS15_glycolysis_gradient_cell_source.csv")
long_source_csv <- file.path(outdir, "FigS15_glycolysis_gradient_long_source.csv")
bin_means_csv <- file.path(outdir, "FigS15_glycolysis_gradient_bin_means_source.csv")
bin_ranges_csv <- file.path(outdir, "FigS15_glycolysis_gradient_bin_ranges_source.csv")
manifest_csv <- file.path(outdir, "FigS15_glycolysis_gradient_manifest.csv")

readr::write_csv(source_wide, cell_source_csv)
readr::write_csv(source_long, long_source_csv)
readr::write_csv(bin_means, bin_means_csv)
readr::write_csv(bin_ranges, bin_ranges_csv)
readr::write_csv(manifest, manifest_csv)

cat("\nSource CSVs written:\n")
cat("  ", cell_source_csv, "\n", sep = "")
cat("  ", long_source_csv, "\n", sep = "")
cat("  ", bin_means_csv, "\n", sep = "")
cat("  ", bin_ranges_csv, "\n", sep = "")
cat("  ", manifest_csv, "\n", sep = "")

# ----------------------------
# Mandatory QC before plotting
# ----------------------------
mean_by_gene_bin <- bin_means %>%
  dplyr::select(gene, glyco_bin_index, mean_log_normalized_expression) %>%
  tidyr::pivot_wider(
    names_from = glyco_bin_index,
    values_from = mean_log_normalized_expression,
    names_prefix = "bin_"
  )

get_bin_mean <- function(gene, bin_index) {
  mean_by_gene_bin %>%
    dplyr::filter(.data$gene == gene) %>%
    dplyr::pull(paste0("bin_", bin_index))
}

high_half_mean <- bin_means %>%
  dplyr::filter(glyco_bin_index >= 6) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(high_half_mean = mean(mean_log_normalized_expression), .groups = "drop")

low_half_mean <- bin_means %>%
  dplyr::filter(glyco_bin_index <= 5) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(low_half_mean = mean(mean_log_normalized_expression), .groups = "drop")

gradient_direction <- high_half_mean %>%
  dplyr::left_join(low_half_mean, by = "gene") %>%
  dplyr::mutate(high_minus_low_half = high_half_mean - low_half_mean)

direction_pass <- all(gradient_direction$high_minus_low_half > 0)

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
    "Target genes present",
    paste(target_genes, collapse = ";"),
    paste(target_genes, collapse = ";"),
    length(missing_target_genes) == 0
  ),
  qc_row(
    "Binning method",
    "10 equal-width bins using cut() over Glycolysis_AUC",
    "10 equal-width AUCell-score bins",
    TRUE
  ),
  qc_row(
    "Number of non-empty bins",
    n_distinct(tumor_hep$glyco_bin_index),
    10,
    n_distinct(tumor_hep$glyco_bin_index) == 10
  ),
  qc_row(
    "No empty equal-width bins",
    paste(bin_counts$n_cells, collapse = ";"),
    "All 10 bins have n_cells > 0",
    all(bin_counts$n_cells > 0)
  ),
  qc_row(
    "Finite bin mean expression values",
    all(is.finite(bin_means$mean_log_normalized_expression)),
    TRUE,
    all(is.finite(bin_means$mean_log_normalized_expression))
  ),
  qc_row(
    "High-half mean expression exceeds low-half mean for all plotted genes",
    paste0(
      gradient_direction$gene,
      "=",
      sprintf("%.4f", gradient_direction$high_minus_low_half),
      collapse = ";"
    ),
    "Positive high-half minus low-half mean for ENO1, LDHA, SPP1, MIF, PTGES",
    direction_pass
  ),
  qc_row(
    "Figure has no manuscript-number title",
    "No ggtitle/labs(title) used",
    "No 'Supplementary Figure S15' title inside figure",
    TRUE
  )
)

qc_status <- ifelse(all(qc$pass), "PASS", "FAIL")
qc <- qc %>% dplyr::mutate(overall_qc = qc_status)

qc_csv <- file.path(outdir, "FigS15_glycolysis_gradient_QC_check.csv")
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
plot_df <- bin_means %>%
  dplyr::mutate(
    gene = factor(gene, levels = target_genes),
    gene_class = factor(
      gene_class,
      levels = c("Glycolytic enzyme", "Immunosuppressive ligand")
    )
  )

fig <- ggplot2::ggplot(
  plot_df,
  ggplot2::aes(
    x = glyco_bin_index,
    y = mean_log_normalized_expression,
    group = gene,
    color = gene,
    linetype = gene_class
  )
) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2.0) +
  ggplot2::scale_x_continuous(
    breaks = 1:10,
    labels = paste0("Bin ", 1:10),
    limits = c(1, 10)
  ) +
  ggplot2::scale_linetype_manual(
    values = c(
      "Glycolytic enzyme" = "solid",
      "Immunosuppressive ligand" = "dashed"
    )
  ) +
  ggplot2::labs(
    x = "Glycolysis activity bin (low to high)",
    y = "Mean log-normalized expression",
    color = "Gene",
    linetype = "Gene class"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
    axis.title = ggplot2::element_text(color = "black"),
    legend.position = "right",
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

pdf_file <- file.path(outdir, "FigS15_glycolysis_gradient.pdf")
png_file <- file.path(outdir, "FigS15_glycolysis_gradient.png")

ggplot2::ggsave(pdf_file, fig, width = 7.2, height = 4.8, useDingbats = FALSE)
ggplot2::ggsave(png_file, fig, width = 7.2, height = 4.8, dpi = 600)

cat("\nQC PASS. Figure files written:\n")
cat("  ", pdf_file, "\n", sep = "")
cat("  ", png_file, "\n", sep = "")

cat("\nSupplementary Figure S15 glycolysis-gradient analysis complete.\n")
