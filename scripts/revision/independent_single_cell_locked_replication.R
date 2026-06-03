#!/usr/bin/env Rscript

# Locked independent single-cell replication workflow.
#
# Intended datasets:
#   GSE189903 (preferred) or GSE146115 (alternative), prepared as a Seurat RDS.
#
# Required environment variables:
#   INDEPENDENT_SEURAT_RDS
#
# Optional environment variables:
#   PATIENT_COLUMN       default: patient
#   SITE_COLUMN          default: site
#   CELLTYPE_COLUMN      default: cell_type
#   TUMOR_LABEL          default: Tumor
#   HEPATOCYTE_LABEL     default: Hepatocyte
#   MALIGNANT_FLAG_COLUMN
#
# The workflow uses the locked 22-gene glycolysis set, fixed four-gene readout
# coefficients, and rank-balanced median grouping. No genes, coefficients, or
# thresholds are optimized in the independent cohort.

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

rds_path <- Sys.getenv("INDEPENDENT_SEURAT_RDS", unset = "independent_hcc_seurat.rds")
patient_col <- Sys.getenv("PATIENT_COLUMN", unset = "patient")
site_col <- Sys.getenv("SITE_COLUMN", unset = "site")
celltype_col <- Sys.getenv("CELLTYPE_COLUMN", unset = "cell_type")
tumor_label <- Sys.getenv("TUMOR_LABEL", unset = "Tumor")
hepatocyte_label <- Sys.getenv("HEPATOCYTE_LABEL", unset = "Hepatocyte")
malignant_flag_col <- Sys.getenv("MALIGNANT_FLAG_COLUMN", unset = "")
outdir <- Sys.getenv(
  "INDEPENDENT_REPLICATION_OUTDIR",
  unset = file.path("results", "independent_single_cell_locked_replication")
)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(rds_path)) stop("Missing independent Seurat object: ", rds_path, call. = FALSE)

glycolysis_22 <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKM", "ALDOA", "ALDOB", "ALDOC",
  "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
  "SLC2A1", "SLC2A3", "PFKFB3", "GCK"
)
four_gene_coef <- c(TPI1 = 0.3041908, ENO1 = 0.9639654, LDHA = 1.3404374, SLC2A1 = 0.2424239)
candidate_ligands <- c("MIF", "SPP1")

get_assay_matrix <- function(seu, assay = "RNA", layer = "data", slot = "data") {
  tryCatch(
    SeuratObject::LayerData(seu, assay = assay, layer = layer),
    error = function(e) Seurat::GetAssayData(seu, assay = assay, slot = slot)
  )
}

rank_balanced_group <- function(df, score_col) {
  df %>%
    arrange(.data[[score_col]], cell) %>%
    mutate(
      score_rank = row_number(),
      GlycoGroup = if_else(score_rank <= ceiling(n() / 2), "GlycoLow", "GlycoHigh")
    )
}

seu <- readRDS(rds_path)
meta <- seu@meta.data %>% rownames_to_column("cell")
required_meta <- c(patient_col, site_col, celltype_col)
missing_meta <- setdiff(required_meta, names(meta))
if (length(missing_meta) > 0) {
  stop("Missing metadata columns: ", paste(missing_meta, collapse = ", "), call. = FALSE)
}

selected <- meta %>%
  filter(
    .data[[site_col]] == tumor_label,
    .data[[celltype_col]] == hepatocyte_label
  ) %>%
  transmute(
    cell,
    patient = as.character(.data[[patient_col]]),
    site = as.character(.data[[site_col]]),
    cell_type = as.character(.data[[celltype_col]])
  )

selection_label <- "tumor-derived hepatocyte-lineage cells"
if (nzchar(malignant_flag_col)) {
  if (!malignant_flag_col %in% names(meta)) {
    stop("MALIGNANT_FLAG_COLUMN not found in metadata: ", malignant_flag_col, call. = FALSE)
  }
  malignant_cells <- meta %>%
    filter(as.logical(.data[[malignant_flag_col]])) %>%
    pull(cell)
  selected <- selected %>% filter(cell %in% malignant_cells)
  selection_label <- "high-confidence malignant hepatocyte-lineage cells"
}

if (nrow(selected) < 100) {
  stop("Too few selected cells for replication: ", nrow(selected), call. = FALSE)
}
if (n_distinct(selected$patient) < 3) {
  stop("Fewer than three independent patients after selection.", call. = FALSE)
}

assay_use <- if ("RNA" %in% names(seu@assays)) "RNA" else DefaultAssay(seu)
expr <- get_assay_matrix(seu, assay = assay_use)
required_genes <- unique(c(glycolysis_22, names(four_gene_coef), candidate_ligands))
missing_genes <- setdiff(required_genes, rownames(expr))
if (length(missing_genes) > 0) {
  stop("Missing required genes: ", paste(missing_genes, collapse = ", "), call. = FALSE)
}

expr_selected <- expr[, selected$cell, drop = FALSE]
rankings <- AUCell_buildRankings(expr_selected, plotStats = FALSE, nCores = 1)
auc <- AUCell_calcAUC(
  list(glycolysis_22 = glycolysis_22),
  rankings,
  aucMaxRank = ceiling(nrow(expr_selected) * 0.05)
)
selected$Glycolysis_AUC <- as.numeric(getAUC(auc)["glycolysis_22", selected$cell])
selected <- rank_balanced_group(selected, "Glycolysis_AUC")

expr_df <- as.matrix(expr_selected[unique(c(names(four_gene_coef), candidate_ligands)), , drop = FALSE]) %>%
  t() %>%
  as.data.frame(check.names = FALSE) %>%
  rownames_to_column("cell") %>%
  mutate(
    four_gene_readout =
      four_gene_coef["TPI1"] * TPI1 +
      four_gene_coef["ENO1"] * ENO1 +
      four_gene_coef["LDHA"] * LDHA +
      four_gene_coef["SLC2A1"] * SLC2A1
  )

analysis_df <- selected %>% left_join(expr_df, by = "cell")
features <- c("MIF", "SPP1", "four_gene_readout")

patient_correlations <- analysis_df %>%
  pivot_longer(all_of(features), names_to = "feature", values_to = "expression") %>%
  group_by(patient, feature) %>%
  summarise(
    n_cells = n(),
    spearman_rho = suppressWarnings(cor(Glycolysis_AUC, expression, method = "spearman")),
    .groups = "drop"
  )

patient_group_means <- analysis_df %>%
  pivot_longer(all_of(features), names_to = "feature", values_to = "expression") %>%
  group_by(patient, GlycoGroup, feature) %>%
  summarise(mean_expression = mean(expression, na.rm = TRUE), .groups = "drop")

patient_effects <- patient_group_means %>%
  pivot_wider(names_from = GlycoGroup, values_from = mean_expression) %>%
  mutate(effect_High_minus_Low = GlycoHigh - GlycoLow)

replication_summary <- patient_effects %>%
  group_by(feature) %>%
  summarise(
    n_patients = n(),
    n_positive_effect = sum(effect_High_minus_Low > 0, na.rm = TRUE),
    mean_effect_High_minus_Low = mean(effect_High_minus_Low, na.rm = TRUE),
    median_effect_High_minus_Low = median(effect_High_minus_Low, na.rm = TRUE),
    paired_wilcox_p = wilcox.test(GlycoHigh, GlycoLow, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(paired_wilcox_FDR = p.adjust(paired_wilcox_p, method = "BH"))

cohort_summary <- tibble(
  selection = selection_label,
  n_cells = nrow(analysis_df),
  n_patients = n_distinct(analysis_df$patient),
  GlycoHigh_n = sum(analysis_df$GlycoGroup == "GlycoHigh"),
  GlycoLow_n = sum(analysis_df$GlycoGroup == "GlycoLow"),
  grouping_rule = "locked rank-balanced median strategy; no re-optimization"
)

write_csv(cohort_summary, file.path(outdir, "independent_replication_cohort_summary.csv"))
write_csv(analysis_df, file.path(outdir, "independent_replication_cell_source.csv"))
write_csv(patient_correlations, file.path(outdir, "independent_replication_patient_correlations.csv"))
write_csv(patient_effects, file.path(outdir, "independent_replication_patient_effects.csv"))
write_csv(replication_summary, file.path(outdir, "independent_replication_summary.csv"))
message("Independent single-cell replication outputs written to: ", outdir)
