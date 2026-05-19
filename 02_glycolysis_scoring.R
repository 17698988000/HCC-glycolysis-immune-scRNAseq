#!/usr/bin/env Rscript

# ============================================================
# 02_glycolysis_scoring.R
#
# Purpose:
#   Define GlycoHigh / GlycoLow tumor-derived hepatocytes using
#   the finalized current-object data convention.
#
# Final manuscript-aligned convention:
#   Input object: seurat_final.rds
#   patient_col  = patient
#   celltype_col = cell_type
#   site_col     = site
#   gly_col      = Glycolysis_AUC
#
# Tumor-derived hepatocytes:
#   site == "Tumor" & cell_type == "Hepatocyte"
#
# GlycoHigh / GlycoLow definition:
#   median split within tumor-derived hepatocytes only
#   High = Glycolysis_AUC > median_cut
#   Low  = all remaining cells
#
# Locked expected values:
#   tumor-derived hepatocytes n = 15,391
#   median Glycolysis_AUC cutoff = 0.2203849
#   GlycoLow  n = 7,696
#   GlycoHigh n = 7,695
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

# ---------- avoid masking ----------
select      <- dplyr::select
mutate      <- dplyr::mutate
filter      <- dplyr::filter
arrange     <- dplyr::arrange
summarise   <- dplyr::summarise
group_by    <- dplyr::group_by
ungroup     <- dplyr::ungroup
count       <- dplyr::count
n_distinct  <- dplyr::n_distinct

# ============================================================
# User-facing inputs
# ============================================================

rds_path <- "seurat_final.rds"

patient_col  <- "patient"
celltype_col <- "cell_type"
site_col     <- "site"
gly_col      <- "Glycolysis_AUC"

outdir <- "results/02_glycolysis_scoring"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

expected_tumor_hep_n <- 15391L
expected_low_n       <- 7696L
expected_high_n      <- 7695L
expected_median_cut  <- 0.2203849
median_tolerance     <- 1e-7

# ============================================================
# Load object
# ============================================================

if (!file.exists(rds_path)) {
  stop("Input file not found: ", rds_path)
}

seu <- readRDS(rds_path)

meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell")

required_cols <- c(patient_col, celltype_col, site_col, gly_col)
missing_cols <- setdiff(required_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    "Required metadata column(s) missing: ",
    paste(missing_cols, collapse = ", "),
    "\nCurrent metadata columns are:\n",
    paste(colnames(meta), collapse = ", ")
  )
}

# Defensive checks against deprecated manuscript-inconsistent columns.
if ("glycolysis_score" %in% colnames(meta)) {
  message(
    "Note: deprecated column 'glycolysis_score' exists, but this script uses the finalized column '",
    gly_col,
    "'."
  )
}

if ("tissue" %in% colnames(meta)) {
  message(
    "Note: deprecated/alternative column 'tissue' exists, but this script uses the finalized column '",
    site_col,
    "'."
  )
}

# ============================================================
# Inspect labels
# ============================================================

site_labels <- sort(unique(as.character(meta[[site_col]])))
celltype_labels <- sort(unique(as.character(meta[[celltype_col]])))

message("Unique site labels:")
message(paste(site_labels, collapse = ", "))

message("Unique cell_type labels:")
message(paste(celltype_labels, collapse = ", "))

if (!"Tumor" %in% site_labels) {
  stop("Expected site label 'Tumor' not found in column: ", site_col)
}

if (!"Hepatocyte" %in% celltype_labels) {
  stop("Expected cell_type label 'Hepatocyte' not found in column: ", celltype_col)
}

if ("Hepatocytes" %in% celltype_labels) {
  message(
    "Note: old label 'Hepatocytes' is present, but the finalized selection uses only 'Hepatocyte'."
  )
}

# ============================================================
# Select tumor-derived hepatocytes
# ============================================================

tumor_hep <- meta %>%
  filter(
    .data[[site_col]] == "Tumor",
    .data[[celltype_col]] == "Hepatocyte"
  ) %>%
  mutate(
    patient = .data[[patient_col]],
    cell_type = .data[[celltype_col]],
    site = .data[[site_col]],
    Glycolysis_AUC = as.numeric(.data[[gly_col]])
  )

if (anyNA(tumor_hep$Glycolysis_AUC)) {
  stop("NA values detected in Glycolysis_AUC among selected tumor-derived hepatocytes.")
}

tumor_hep_n <- nrow(tumor_hep)

message("Selected tumor-derived hepatocytes: ", tumor_hep_n)

if (tumor_hep_n != expected_tumor_hep_n) {
  stop(
    "Selected tumor-derived hepatocytes != ",
    expected_tumor_hep_n,
    ". Observed n = ",
    tumor_hep_n,
    ". Check site/cell_type labels and input object."
  )
}

# ============================================================
# Median split within tumor-derived hepatocytes
# ============================================================

median_cut <- median(tumor_hep$Glycolysis_AUC, na.rm = TRUE)

message("Median Glycolysis_AUC cutoff: ", format(median_cut, digits = 10))

if (abs(median_cut - expected_median_cut) > median_tolerance) {
  stop(
    "Median Glycolysis_AUC cutoff does not match locked value. ",
    "Expected ",
    expected_median_cut,
    "; observed ",
    format(median_cut, digits = 10),
    "."
  )
}

tumor_hep <- tumor_hep %>%
  mutate(
    gly_group = ifelse(Glycolysis_AUC > median_cut, "GlycoHigh", "GlycoLow"),
    gly_group = factor(gly_group, levels = c("GlycoLow", "GlycoHigh"))
  )

group_counts <- table(tumor_hep$gly_group)

low_n <- as.integer(group_counts[["GlycoLow"]])
high_n <- as.integer(group_counts[["GlycoHigh"]])

message("GlycoLow n: ", low_n)
message("GlycoHigh n: ", high_n)

if (low_n != expected_low_n || high_n != expected_high_n) {
  stop(
    "Glyco group counts do not match locked values. ",
    "Expected GlycoLow = ",
    expected_low_n,
    ", GlycoHigh = ",
    expected_high_n,
    "; observed GlycoLow = ",
    low_n,
    ", GlycoHigh = ",
    high_n,
    "."
  )
}

# ============================================================
# Write source outputs
# ============================================================

glyco_assignment <- tumor_hep %>%
  select(
    cell,
    patient,
    site,
    cell_type,
    Glycolysis_AUC,
    gly_group
  ) %>%
  arrange(patient, gly_group, cell)

patient_summary <- glyco_assignment %>%
  group_by(patient, gly_group) %>%
  summarise(
    n_cells = dplyr::n(),
    mean_Glycolysis_AUC = mean(Glycolysis_AUC, na.rm = TRUE),
    median_Glycolysis_AUC = median(Glycolysis_AUC, na.rm = TRUE),
    min_Glycolysis_AUC = min(Glycolysis_AUC, na.rm = TRUE),
    max_Glycolysis_AUC = max(Glycolysis_AUC, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(patient, gly_group)

qc_check <- tibble::tibble(
  check_name = c(
    "input_rds",
    "patient_col",
    "celltype_col",
    "site_col",
    "gly_col",
    "tumor_hepatocyte_rule",
    "tumor_hepatocyte_n",
    "median_Glycolysis_AUC_cutoff",
    "GlycoLow_n",
    "GlycoHigh_n",
    "High_definition",
    "Low_definition",
    "overall_status"
  ),
  expected = c(
    rds_path,
    "patient",
    "cell_type",
    "site",
    "Glycolysis_AUC",
    'site == "Tumor" & cell_type == "Hepatocyte"',
    as.character(expected_tumor_hep_n),
    as.character(expected_median_cut),
    as.character(expected_low_n),
    as.character(expected_high_n),
    "Glycolysis_AUC > median_cut",
    "Glycolysis_AUC <= median_cut",
    "PASS"
  ),
  observed = c(
    rds_path,
    patient_col,
    celltype_col,
    site_col,
    gly_col,
    'site == "Tumor" & cell_type == "Hepatocyte"',
    as.character(tumor_hep_n),
    format(median_cut, digits = 10),
    as.character(low_n),
    as.character(high_n),
    "Glycolysis_AUC > median_cut",
    "Glycolysis_AUC <= median_cut",
    "PASS"
  ),
  pass = c(
    file.exists(rds_path),
    patient_col == "patient",
    celltype_col == "cell_type",
    site_col == "site",
    gly_col == "Glycolysis_AUC",
    TRUE,
    tumor_hep_n == expected_tumor_hep_n,
    abs(median_cut - expected_median_cut) <= median_tolerance,
    low_n == expected_low_n,
    high_n == expected_high_n,
    TRUE,
    TRUE,
    TRUE
  )
)

write_csv(
  glyco_assignment,
  file.path(outdir, "02_glycolysis_scoring_tumor_hepatocyte_glyco_assignment.csv")
)

write_csv(
  patient_summary,
  file.path(outdir, "02_glycolysis_scoring_patient_summary.csv")
)

write_csv(
  qc_check,
  file.path(outdir, "02_glycolysis_scoring_QC_check.csv")
)

# Optional: save a compact RDS assignment object for downstream scripts.
saveRDS(
  glyco_assignment,
  file.path(outdir, "02_glycolysis_scoring_tumor_hepatocyte_glyco_assignment.rds")
)

# Optional: add glyco group metadata back to the Seurat object without changing
# non-selected cells. This is useful for downstream scripts that expect metadata.
seu$TumorHep_GlycoGroup <- NA_character_
seu$TumorHep_GlycoGroup[match(glyco_assignment$cell, colnames(seu))] <-
  as.character(glyco_assignment$gly_group)

seu$TumorHep_Glycolysis_AUC_MedianCut <- NA_real_
seu$TumorHep_Glycolysis_AUC_MedianCut[match(glyco_assignment$cell, colnames(seu))] <-
  median_cut

saveRDS(
  seu,
  file.path(outdir, "seurat_final_with_TumorHep_GlycoGroup.rds")
)

message("QC PASS.")
message("Outputs written to: ", outdir)
message("Primary assignment CSV: 02_glycolysis_scoring_tumor_hepatocyte_glyco_assignment.csv")
message("QC CSV: 02_glycolysis_scoring_QC_check.csv")
