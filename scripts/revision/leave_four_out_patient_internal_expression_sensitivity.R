#!/usr/bin/env Rscript

# Leave-four-out and patient-internal grouping sensitivity analysis.
#
# Required input:
#   SEURAT_FINAL_RDS: Seurat object containing patient, site, cell_type metadata.
#
# This script removes TPI1, ENO1, LDHA, and SLC2A1 from the locked 22-gene
# glycolysis set, recalculates AUCell scores, applies the locked rank-balanced
# median-split strategy, and evaluates patient-level expression effects for the
# four-gene readout and candidate ligands. It does not re-optimize genes or
# thresholds.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(AUCell)
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

set.seed(20260603)

rds_path <- Sys.getenv("SEURAT_FINAL_RDS", unset = "seurat_final.rds")
outdir <- Sys.getenv(
  "REVISION_SENSITIVITY_OUTDIR",
  unset = file.path("results", "leave_four_out_patient_internal_sensitivity")
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(rds_path)) {
  stop("Missing Seurat object: ", rds_path, call. = FALSE)
}

glycolysis_22 <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKM", "ALDOA", "ALDOB", "ALDOC",
  "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
  "SLC2A1", "SLC2A3", "PFKFB3", "GCK"
)
readout_genes <- c("TPI1", "ENO1", "LDHA", "SLC2A1")
leave_four_out_genes <- setdiff(glycolysis_22, readout_genes)
association_genes <- c(readout_genes, "MIF", "SPP1", "PTGES")

rank_balanced_group <- function(df, score_col, patient_internal = FALSE) {
  if (patient_internal) {
    df %>%
      group_by(patient) %>%
      arrange(.data[[score_col]], cell, .by_group = TRUE) %>%
      mutate(
        rank_in_group = row_number(),
        n_in_group = n(),
        group = if_else(rank_in_group <= ceiling(n_in_group / 2), "GlycoLow", "GlycoHigh")
      ) %>%
      ungroup()
  } else {
    df %>%
      arrange(.data[[score_col]], cell) %>%
      mutate(
        rank_in_group = row_number(),
        n_in_group = n(),
        group = if_else(rank_in_group <= ceiling(n_in_group / 2), "GlycoLow", "GlycoHigh")
      )
  }
}

get_assay_matrix <- function(seu, assay = "RNA", layer = "data", slot = "data") {
  tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = slot)
  )
}

seu <- readRDS(rds_path)
meta <- seu@meta.data %>% rownames_to_column("cell")
required_meta <- c("patient", "site", "cell_type")
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

tumor_hep <- meta %>%
  filter(site == "Tumor", cell_type == "Hepatocyte") %>%
  select(cell, patient, site, cell_type)
if (nrow(tumor_hep) != 15391) {
  stop("Expected 15,391 tumor-derived hepatocytes; observed ", nrow(tumor_hep), call. = FALSE)
}

assay_use <- if ("RNA" %in% names(seu@assays)) "RNA" else DefaultAssay(seu)
expr <- get_assay_matrix(seu, assay = assay_use)
missing_score_genes <- setdiff(leave_four_out_genes, rownames(expr))
if (length(missing_score_genes) > 0) {
  stop("Missing leave-four-out score genes: ", paste(missing_score_genes, collapse = ", "), call. = FALSE)
}
missing_assoc_genes <- setdiff(association_genes, rownames(expr))
if (length(missing_assoc_genes) > 0) {
  stop("Missing association genes: ", paste(missing_assoc_genes, collapse = ", "), call. = FALSE)
}

expr_tumor_hep <- expr[, tumor_hep$cell, drop = FALSE]
rankings <- AUCell_buildRankings(expr_tumor_hep, plotStats = FALSE, nCores = 1)
auc <- AUCell_calcAUC(
  list(glycolysis_leave_four_out = leave_four_out_genes),
  rankings,
  aucMaxRank = ceiling(nrow(expr_tumor_hep) * 0.05)
)
score <- as.numeric(getAUC(auc)["glycolysis_leave_four_out", tumor_hep$cell])
names(score) <- tumor_hep$cell
tumor_hep$leave_four_out_AUC <- score[tumor_hep$cell]

global_assignment <- rank_balanced_group(tumor_hep, "leave_four_out_AUC", FALSE) %>%
  select(cell, patient, leave_four_out_AUC, leave_four_out_global_group = group)
patient_assignment <- rank_balanced_group(tumor_hep, "leave_four_out_AUC", TRUE) %>%
  select(cell, patient, leave_four_out_AUC, leave_four_out_patient_internal_group = group)
assignments <- global_assignment %>%
  left_join(
    patient_assignment %>% select(-patient, -leave_four_out_AUC),
    by = "cell"
  )

expr_df <- as.matrix(expr_tumor_hep[association_genes, , drop = FALSE]) %>%
  t() %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("cell")

analyze_grouping <- function(group_col, grouping_name) {
  joined <- assignments %>%
    select(cell, patient, group = all_of(group_col)) %>%
    left_join(expr_df, by = "cell")
  means <- joined %>%
    pivot_longer(all_of(association_genes), names_to = "gene", values_to = "expression") %>%
    group_by(patient, group, gene) %>%
    summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop")
  effects <- means %>%
    pivot_wider(names_from = group, values_from = mean_expression) %>%
    mutate(effect_High_minus_Low = GlycoHigh - GlycoLow, grouping = grouping_name)
  summary <- effects %>%
    group_by(grouping, gene) %>%
    summarise(
      n_patients = n(),
      n_positive = sum(effect_High_minus_Low > 0, na.rm = TRUE),
      mean_effect_High_minus_Low = mean(effect_High_minus_Low, na.rm = TRUE),
      median_effect_High_minus_Low = median(effect_High_minus_Low, na.rm = TRUE),
      paired_wilcox_p = wilcox.test(GlycoHigh, GlycoLow, paired = TRUE, exact = FALSE)$p.value,
      .groups = "drop"
    ) %>%
    mutate(paired_wilcox_FDR = p.adjust(paired_wilcox_p, method = "BH"))
  list(means = means, effects = effects, summary = summary)
}

global_results <- analyze_grouping("leave_four_out_global_group", "leave_four_out_global")
patient_results <- analyze_grouping(
  "leave_four_out_patient_internal_group",
  "leave_four_out_patient_internal"
)

continuous_correlations <- assignments %>%
  select(cell, patient, leave_four_out_AUC) %>%
  left_join(expr_df, by = "cell") %>%
  pivot_longer(all_of(association_genes), names_to = "gene", values_to = "expression") %>%
  group_by(patient, gene) %>%
  summarise(
    spearman_rho = suppressWarnings(cor(leave_four_out_AUC, expression, method = "spearman")),
    .groups = "drop"
  )

write_csv(assignments, file.path(outdir, "leave_four_out_assignments.csv"))
write_csv(bind_rows(global_results$means, patient_results$means), file.path(outdir, "leave_four_out_patient_means.csv"))
write_csv(bind_rows(global_results$effects, patient_results$effects), file.path(outdir, "leave_four_out_patient_effects.csv"))
write_csv(bind_rows(global_results$summary, patient_results$summary), file.path(outdir, "leave_four_out_patient_summary.csv"))
write_csv(continuous_correlations, file.path(outdir, "leave_four_out_continuous_correlations.csv"))

message("Leave-four-out sensitivity outputs written to: ", outdir)
