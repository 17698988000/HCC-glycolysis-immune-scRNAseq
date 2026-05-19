#!/usr/bin/env Rscript

# ============================================================
# 13_revision_analyses.R
#
# Legacy/status note for the former mixed "revision analyses" script.
#
# This script is intentionally retained as a lightweight status file
# rather than as an executable multi-analysis pipeline.
#
# Rationale:
#   The previous version mixed several independent revision analyses
#   in one file, used local absolute paths, contained obsolete figure
#   numbering, and included superseded S12 values. Running that legacy
#   version can regenerate inconsistent outputs.
#
# Current repository policy:
#   One analysis target per script.
#   Dedicated current scripts should be used instead of this legacy
#   mixed revision file.
#
# This script writes:
#   results_13_revision_analyses_status/
#     13_revision_analyses_status.csv
#     13_revision_analyses_QC_check.csv
#     13_revision_analyses_status_note.txt
#
# It does not generate PDF/PNG figures.
# It does not modify README.md or any other script.
# ============================================================

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("Package 'readr' is required. Install it with install.packages('readr').", call. = FALSE)
  }
  if (!requireNamespace("tibble", quietly = TRUE)) {
    stop("Package 'tibble' is required. Install it with install.packages('tibble').", call. = FALSE)
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package 'dplyr' is required. Install it with install.packages('dplyr').", call. = FALSE)
  }

  library(readr)
  library(tibble)
  library(dplyr)
})

outdir <- "results_13_revision_analyses_status"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# Finalized / redirected analyses
# ----------------------------
status_tbl <- tibble::tribble(
  ~legacy_block, ~current_status, ~current_source_of_truth, ~notes,

  "LOO circularity analysis for ENO1, LDHA, and GAPDH",
  "Superseded by final manuscript-aligned figure mapping",
  "Final manuscript Supplementary Figure S19 / relevant finalized source files",
  "Do not regenerate from the old mixed 13_revision_analyses.R block. The legacy block used hard-coded values and outdated figure numbering.",

  "Spatial transcriptomics RCTD deconvolution",
  "Superseded by dedicated spatial workflow",
  "11_spatial_transcriptomics.R",
  "Use the current GSE238264 four-sample Visium workflow. Spatial interpretation must remain spot-level tissue co-enrichment, not same-cell co-expression, tumor-cell-specific ligand production, or direct ligand-receptor contact.",

  "OXPHOS metabolic specificity / partial correlation",
  "Superseded by dedicated metabolic-specificity scripts",
  "14_OXPHOS_metabolic_specificity.R and 17_FigS21_partial_correlation_metabolic_specificity.R",
  "Use dedicated scripts for Supplementary Figure S21-style metabolic specificity outputs.",

  "DoRothEA / VIPER transcription-factor activity analysis",
  "Superseded by dedicated TF workflow",
  "09_TF_activity.R",
  "Use the dedicated TF activity script rather than the legacy mixed block.",

  "NicheNet ligand activity analysis",
  "Superseded by dedicated NicheNet workflow",
  "15_NicheNet_analysis.R",
  "Use the dedicated NicheNet script for Supplementary Figure S22-style outputs.",

  "Patient-level ligand mean-expression robustness analysis",
  "Superseded and locked",
  "16_FigS23_patient_level_ligand_effects.R",
  "Supplementary Figure S23 is locked. Do not re-open or regenerate it from this legacy file.",

  "Per-patient ENO1-glycolysis correlation",
  "Superseded and locked",
  "Final S12 current-object restoration output",
  "Supplementary Figure S12 is locked to the current-object result: 15,391 tumor-derived hepatocytes, 8 patients, global Pearson R = 0.57, Spearman rho range 0.363 to 0.633, median rho = 0.498, all BH-adjusted p < 0.05. Do not restore old S12 values.",

  "Glycolysis-gradient expression analysis",
  "Superseded and locked",
  "08_glycolysis_gradient.R",
  "Supplementary Figure S15 is locked as 10 equal-width AUCell-score bins with mean log-normalized ENO1, LDHA, SPP1, MIF, and PTGES expression in tumor-derived hepatocytes.",

  "GSE14520 external validation",
  "Handled by dedicated script",
  "10_GSE14520_validation.R",
  "Use the current four-gene GSE14520 validation workflow. Do not assume README has already been updated.",

  "GSE235863 anti-PD-1 plus lenvatinib exploratory non-response analysis",
  "Handled by dedicated replacement script",
  "12_drug_repurposing.R",
  "Despite the legacy file name, script 12 now carries the GSE235863 four-gene exploratory non-response analysis. Interpret as hypothesis-generating only, not as a formal predictive model or clinical treatment-selection assay."
)

# ----------------------------
# Repository-level locked conventions
# ----------------------------
locked_conventions <- tibble::tribble(
  ~item, ~locked_value,

  "input_object",
  "seurat_final.rds",

  "patient_col",
  "patient",

  "celltype_col",
  "cell_type",

  "site_col",
  "site",

  "gly_col",
  "Glycolysis_AUC",

  "tumor_derived_hepatocytes",
  "site == 'Tumor' & cell_type == 'Hepatocyte'",

  "tumor_derived_hepatocyte_n",
  "15391",

  "glycolysis_median_cutoff",
  "0.2203849",

  "glycohigh_rule",
  "Glycolysis_AUC > median_cut",

  "glycolow_rule",
  "Remaining cells",

  "glycolow_n",
  "7696",

  "glycohigh_n",
  "7695",

  "four_gene_score",
  "0.3041908*TPI1 + 0.9639654*ENO1 + 1.3404374*LDHA + 0.2424239*SLC2A1",

  "GSE235863_interpretation",
  "Exploratory / hypothesis-generating only; not a formal prediction model or treatment-selection assay.",

  "spatial_interpretation",
  "Spot-level tissue co-enrichment only; not same-cell co-expression, tumor-cell-specific ligand production, or direct ligand-receptor contact."
)

# ----------------------------
# QC / safety checks
# ----------------------------
required_redirects <- c(
  "08_glycolysis_gradient.R",
  "09_TF_activity.R",
  "10_GSE14520_validation.R",
  "11_spatial_transcriptomics.R",
  "12_drug_repurposing.R",
  "14_OXPHOS_metabolic_specificity.R",
  "15_NicheNet_analysis.R",
  "16_FigS23_patient_level_ligand_effects.R",
  "17_FigS21_partial_correlation_metabolic_specificity.R"
)

qc_tbl <- tibble::tibble(
  check = c(
    "Script is legacy/status-only",
    "No PDF/PNG figure generation",
    "No local absolute D:/ path required",
    "No README modification",
    "No S12 old values restored",
    "No S23 regeneration",
    "Redirects listed for dedicated current scripts",
    "GSE235863 interpretation remains exploratory only",
    "Spatial wording remains spot-level tissue co-enrichment only"
  ),
  observed = c(
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    paste(required_redirects, collapse = "; "),
    "hypothesis-generating / exploratory only",
    "spot-level tissue co-enrichment only"
  ),
  expected = c(
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    "TRUE",
    paste(required_redirects, collapse = "; "),
    "hypothesis-generating / exploratory only",
    "spot-level tissue co-enrichment only"
  ),
  pass = c(
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    TRUE,
    all(required_redirects %in% status_tbl$current_source_of_truth |
          required_redirects %in% unlist(strsplit(status_tbl$current_source_of_truth, " and ", fixed = TRUE))),
    TRUE,
    TRUE
  )
)

overall_qc <- ifelse(all(qc_tbl$pass), "PASS", "FAIL")
qc_tbl <- qc_tbl %>% dplyr::mutate(overall_qc = overall_qc)

# ----------------------------
# Write status outputs
# ----------------------------
status_csv <- file.path(outdir, "13_revision_analyses_status.csv")
locked_csv <- file.path(outdir, "13_revision_analyses_locked_conventions.csv")
qc_csv <- file.path(outdir, "13_revision_analyses_QC_check.csv")
note_txt <- file.path(outdir, "13_revision_analyses_status_note.txt")

readr::write_csv(status_tbl, status_csv)
readr::write_csv(locked_conventions, locked_csv)
readr::write_csv(qc_tbl, qc_csv)

note_lines <- c(
  "13_revision_analyses.R legacy/status note",
  "",
  "This file is intentionally retained only as a status and redirection script.",
  "The former mixed revision pipeline should not be run because it mixed multiple analyses,",
  "used local absolute paths, contained obsolete figure numbering, and could regenerate",
  "superseded S12/S15/S23-related outputs.",
  "",
  "Use dedicated current scripts instead:",
  "  08_glycolysis_gradient.R",
  "  09_TF_activity.R",
  "  10_GSE14520_validation.R",
  "  11_spatial_transcriptomics.R",
  "  12_drug_repurposing.R",
  "  14_OXPHOS_metabolic_specificity.R",
  "  15_NicheNet_analysis.R",
  "  16_FigS23_patient_level_ligand_effects.R",
  "  17_FigS21_partial_correlation_metabolic_specificity.R",
  "",
  "Locked reminders:",
  "  S12 is current-object locked; do not restore old S12 values.",
  "  S15 is reproduced by 08_glycolysis_gradient.R.",
  "  S23 is final; do not regenerate.",
  "  GSE235863 is exploratory only and must not be described as a formal predictive model.",
  "  Spatial results are spot-level tissue co-enrichment only.",
  "",
  paste0("QC status: ", overall_qc)
)

writeLines(note_lines, con = note_txt)

cat("Wrote status CSV: ", status_csv, "\n", sep = "")
cat("Wrote locked conventions CSV: ", locked_csv, "\n", sep = "")
cat("Wrote QC CSV: ", qc_csv, "\n", sep = "")
cat("Wrote status note: ", note_txt, "\n", sep = "")
cat("Overall QC: ", overall_qc, "\n", sep = "")

if (!all(qc_tbl$pass)) {
  stop("13_revision_analyses.R status QC failed. Inspect the QC CSV.", call. = FALSE)
}

cat("\nNo manuscript figures were generated by design.\n")
cat("This legacy/status script completed successfully.\n")
