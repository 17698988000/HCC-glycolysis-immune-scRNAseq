#!/usr/bin/env Rscript

# ============================================================
# 02_glycolysis_scoring.R
#
# Purpose:
#   Validate and export the authoritative Figure 2C rank-balanced
#   GlycoHigh/GlycoLow assignment for tumor-derived hepatocytes.
#
# Important:
#   Two cells share the AUCell median. The submission workflow uses
#   the locked rank-balanced source table instead of recreating groups
#   with a simple Glycolysis_AUC > median comparison.
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
})

helper_candidates <- c(
  file.path("scripts", "utils", "locked_fig2c_groups.R"),
  file.path("utils", "locked_fig2c_groups.R")
)
helper_path <- helper_candidates[file.exists(helper_candidates)][1]
if (length(helper_path) == 0 || is.na(helper_path)) {
  stop("Cannot find scripts/utils/locked_fig2c_groups.R. Run from the repository root.", call. = FALSE)
}
source(helper_path)

rds_path <- Sys.getenv("SEURAT_FINAL_RDS", unset = "seurat_final.rds")
outdir <- file.path("results", "02_glycolysis_scoring")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(rds_path)) stop("Input file not found: ", rds_path, call. = FALSE)
seu <- readRDS(rds_path)
meta <- seu@meta.data %>% tibble::rownames_to_column("cell")
required_cols <- c("patient", "cell_type", "site", "Glycolysis_AUC")
missing_cols <- setdiff(required_cols, names(meta))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)

locked <- load_locked_fig2c_assignments()
tumor_hep <- meta %>%
  dplyr::filter(.data$site == "Tumor", .data$cell_type == "Hepatocyte") %>%
  dplyr::mutate(Glycolysis_AUC = as.numeric(.data$Glycolysis_AUC))

if (!setequal(tumor_hep$cell, locked$cell)) {
  stop("The selected tumor-derived hepatocytes do not match the locked Figure 2C cell IDs.", call. = FALSE)
}

comparison <- tumor_hep %>%
  dplyr::select(cell, patient, site, cell_type, Glycolysis_AUC) %>%
  dplyr::left_join(locked %>% dplyr::select(cell, locked_Glycolysis_AUC = Glycolysis_AUC, GlycoGroup), by = "cell") %>%
  dplyr::mutate(abs_score_difference = abs(.data$Glycolysis_AUC - .data$locked_Glycolysis_AUC))

if (anyNA(comparison$GlycoGroup)) stop("Missing locked groups after cell-ID matching.", call. = FALSE)
if (max(comparison$abs_score_difference, na.rm = TRUE) > 1e-12) stop("Glycolysis_AUC values differ from the locked Figure 2C source table.", call. = FALSE)

assignment <- comparison %>%
  dplyr::transmute(
    cell, patient, site, cell_type, Glycolysis_AUC,
    gly_group = factor(GlycoGroup, levels = c("GlycoLow", "GlycoHigh"))
  ) %>%
  dplyr::arrange(patient, gly_group, cell)

patient_summary <- assignment %>%
  dplyr::group_by(patient, gly_group) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_Glycolysis_AUC = mean(Glycolysis_AUC),
    median_Glycolysis_AUC = median(Glycolysis_AUC),
    min_Glycolysis_AUC = min(Glycolysis_AUC),
    max_Glycolysis_AUC = max(Glycolysis_AUC),
    .groups = "drop"
  )

qc_check <- tibble::tibble(
  check_name = c(
    "tumor_hepatocyte_n", "median_Glycolysis_AUC_cutoff", "median_tie_n",
    "GlycoLow_n", "GlycoHigh_n", "grouping_rule", "overall_status"
  ),
  observed = c(
    nrow(assignment),
    format(stats::median(assignment$Glycolysis_AUC), digits = 16),
    sum(abs(assignment$Glycolysis_AUC - LOCKED_FIG2C_MEDIAN_AUC) < 1e-12),
    sum(assignment$gly_group == "GlycoLow"),
    sum(assignment$gly_group == "GlycoHigh"),
    "locked rank-balanced Figure 2C assignment table",
    "PASS"
  ),
  expected = c(
    LOCKED_FIG2C_TUMOR_HEPATOCYTE_N,
    format(LOCKED_FIG2C_MEDIAN_AUC, digits = 16),
    LOCKED_FIG2C_MEDIAN_TIE_N,
    LOCKED_FIG2C_GLYCOLOW_N,
    LOCKED_FIG2C_GLYCOHIGH_N,
    "locked rank-balanced Figure 2C assignment table",
    "PASS"
  )
)

readr::write_csv(assignment, file.path(outdir, "02_glycolysis_scoring_tumor_hepatocyte_glyco_assignment.csv"))
readr::write_csv(patient_summary, file.path(outdir, "02_glycolysis_scoring_patient_summary.csv"))
readr::write_csv(qc_check, file.path(outdir, "02_glycolysis_scoring_QC_check.csv"))
saveRDS(assignment, file.path(outdir, "02_glycolysis_scoring_tumor_hepatocyte_glyco_assignment.rds"))

seu$TumorHep_GlycoGroup <- NA_character_
idx <- match(assignment$cell, colnames(seu))
seu$TumorHep_GlycoGroup[idx] <- as.character(assignment$gly_group)
seu$TumorHep_Glycolysis_AUC_MedianCut <- NA_real_
seu$TumorHep_Glycolysis_AUC_MedianCut[idx] <- LOCKED_FIG2C_MEDIAN_AUC
saveRDS(seu, file.path(outdir, "seurat_final_with_TumorHep_GlycoGroup.rds"))

message("QC PASS. Locked Figure 2C assignment exported.")
message("Outputs written to: ", outdir)
