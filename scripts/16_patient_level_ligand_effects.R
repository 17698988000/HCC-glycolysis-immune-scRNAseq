# ============================================================
# Supplementary Figure S23
# Patient-level GlycoHigh-minus-GlycoLow mean-expression effects
# for SPP1, MIF, and PTGES in tumor-derived hepatocytes
#
# Required object:
#   seurat_final.rds
#
# Fixed columns:
#   patient_col  = "patient"
#   celltype_col = "cell_type"
#   site_col     = "site"
#   gly_col      = "Glycolysis_AUC"
#
# Expected selected tumor-derived hepatocytes:
#   n = 15,391
#
# Expected GlycoHigh/GlycoLow split:
#   Low  = 7,696
#   High = 7,695
#
# Final S23 locked QC targets:
#   MIF mean effect High-Low   ~= 0.388
#   MIF median effect High-Low ~= 0.246
#   MIF paired Wilcoxon p      ~= 0.0300
#   MIF paired Wilcoxon FDR    ~= 0.0432
#   SPP1 paired Wilcoxon FDR   ~= 0.0432
#   PTGES paired Wilcoxon FDR  ~= 0.2084
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
})


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

# -------------------------
# User paths
# -------------------------

# Recommended use:
#   Put seurat_final.rds in the working directory, or set:
#   Sys.setenv(SEURAT_FINAL_RDS = "path/to/seurat_final.rds")
#
# Command line example:
#   Rscript 16_FigS23_patient_level_ligand_effects.R

rds_path <- Sys.getenv("SEURAT_FINAL_RDS", unset = "seurat_final.rds")

outdir <- file.path("results", "FigS23_patient_level_ligand_effects")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# -------------------------
# Fixed metadata columns
# -------------------------

patient_col  <- "patient"
celltype_col <- "cell_type"
site_col     <- "site"
gly_col      <- "Glycolysis_AUC"

# -------------------------
# Load Seurat object
# -------------------------

if (!file.exists(rds_path)) {
  stop(
    "Cannot find seurat_final.rds. Current rds_path = ", rds_path,
    "\nPlace seurat_final.rds in the working directory or set SEURAT_FINAL_RDS."
  )
}

seu <- readRDS(rds_path)
meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell")

required_cols <- c(patient_col, celltype_col, site_col, gly_col)

missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

# -------------------------
# Get normalized expression matrix
# Seurat v5/v4 compatible
# -------------------------

assay_use <- if ("RNA" %in% names(seu@assays)) "RNA" else DefaultAssay(seu)

get_norm_mat <- function(seu, assay = assay_use) {
  tryCatch(
    {
      SeuratObject::LayerData(seu, assay = assay, layer = "data")
    },
    error = function(e) {
      Seurat::GetAssayData(seu, assay = assay, slot = "data")
    }
  )
}

norm_mat <- get_norm_mat(seu, assay_use)

# -------------------------
# Select tumor-derived hepatocytes
# -------------------------

tumor_hepatocyte_labels <- c("Hepatocyte")

tumor_hep <- meta %>%
  dplyr::filter(
    .data[[site_col]] == "Tumor",
    .data[[celltype_col]] %in% tumor_hepatocyte_labels
  ) %>%
  dplyr::mutate(
    patient = .data[[patient_col]],
    Glycolysis_AUC = .data[[gly_col]]
  )

cat("\nSelected tumor-derived hepatocytes:", nrow(tumor_hep), "\n")

if (nrow(tumor_hep) != 15391) {
  cat("\nAvailable Tumor-site cell_type labels:\n")
  print(sort(table(meta[[celltype_col]][meta[[site_col]] == "Tumor"]), decreasing = TRUE))
  stop(
    "Selected tumor-derived hepatocytes != 15,391. ",
    "Please check tumor_hepatocyte_labels and metadata labels."
  )
}

# -------------------------
# GlycoHigh / GlycoLow from authoritative rank-balanced Figure 2C assignment
# -------------------------

median_cut <- LOCKED_FIG2C_MEDIAN_AUC
tumor_hep <- tumor_hep %>%
  dplyr::mutate(
    gly_group = map_locked_fig2c_groups(cell, high_label = "High", low_label = "Low"),
    gly_group = factor(gly_group, levels = c("Low", "High"))
  )
if (anyNA(tumor_hep$gly_group)) stop("Locked Figure 2C grouping could not be mapped to every tumor-derived hepatocyte.")

cat("Median Glycolysis_AUC cutoff:", median_cut, "\n")
cat("\nGlyco group counts:\n")
print(table(tumor_hep$gly_group))

gly_table <- table(tumor_hep$gly_group)

if (!all(c("Low", "High") %in% names(gly_table))) {
  stop("GlycoHigh/GlycoLow split failed.")
}

if (as.integer(gly_table[["Low"]]) != 7696 || as.integer(gly_table[["High"]]) != 7695) {
  stop("Unexpected Glyco group counts. Expected Low = 7696 and High = 7695.")
}

# -------------------------
# Helper: pull normalized gene expression
# -------------------------

pull_gene_expr <- function(genes, cells, norm_mat) {
  genes_present <- intersect(genes, rownames(norm_mat))
  
  if (length(genes_present) != length(genes)) {
    warning("Missing genes: ", paste(setdiff(genes, genes_present), collapse = ", "))
  }
  
  if (length(genes_present) == 0) {
    stop("None of the requested genes were found in the expression matrix.")
  }
  
  expr <- as.matrix(norm_mat[genes_present, cells, drop = FALSE])
  expr <- t(expr)
  expr <- as.data.frame(expr)
  expr <- tibble::rownames_to_column(expr, "cell")
  
  expr
}

# -------------------------
# S23 analysis
# -------------------------

genes_s23 <- c("SPP1", "MIF", "PTGES")

s23_expr <- pull_gene_expr(genes_s23, tumor_hep$cell, norm_mat)

s23_df <- tumor_hep %>%
  dplyr::select(cell, patient, gly_group) %>%
  dplyr::left_join(s23_expr, by = "cell")

# Patient-level means
s23_means <- s23_df %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(genes_s23),
    names_to = "gene",
    values_to = "expression"
  ) %>%
  dplyr::group_by(patient, gly_group, gene) %>%
  dplyr::summarise(
    mean_expression = mean(expression, na.rm = TRUE),
    .groups = "drop"
  )

# High-minus-Low effects
s23_effects <- s23_means %>%
  tidyr::pivot_wider(
    names_from = gly_group,
    values_from = mean_expression
  ) %>%
  dplyr::mutate(
    effect_High_minus_Low = High - Low
  ) %>%
  dplyr::arrange(gene, patient)

# Paired patient-level Wilcoxon tests
s23_stats <- s23_effects %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_patients = dplyr::n(),
    mean_low = mean(Low, na.rm = TRUE),
    mean_high = mean(High, na.rm = TRUE),
    mean_effect_High_minus_Low = mean(effect_High_minus_Low, na.rm = TRUE),
    median_effect_High_minus_Low = median(effect_High_minus_Low, na.rm = TRUE),
    paired_wilcox_p = wilcox.test(High, Low, paired = TRUE, exact = FALSE)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    paired_wilcox_FDR = p.adjust(paired_wilcox_p, method = "BH")
  )

# Save source tables
readr::write_csv(
  s23_means,
  file.path(outdir, "FigS23_patient_level_means.csv")
)

readr::write_csv(
  s23_effects,
  file.path(outdir, "FigS23_patient_level_effects.csv")
)

readr::write_csv(
  s23_stats,
  file.path(outdir, "FigS23_patient_level_stats.csv")
)

cat("\n==================== S23 STATS ====================\n")
print(s23_stats)

# -------------------------
# Automatic QC against locked final text/caption values
# -------------------------

get_stat <- function(g, colname) {
  s23_stats %>%
    dplyr::filter(gene == g) %>%
    dplyr::pull(dplyr::all_of(colname))
}

qc_table <- tibble::tibble(
  item = c(
    "MIF mean effect High-Low",
    "MIF median effect High-Low",
    "MIF paired Wilcoxon p",
    "MIF paired Wilcoxon FDR",
    "SPP1 paired Wilcoxon FDR",
    "PTGES paired Wilcoxon FDR"
  ),
  observed = c(
    get_stat("MIF", "mean_effect_High_minus_Low"),
    get_stat("MIF", "median_effect_High_minus_Low"),
    get_stat("MIF", "paired_wilcox_p"),
    get_stat("MIF", "paired_wilcox_FDR"),
    get_stat("SPP1", "paired_wilcox_FDR"),
    get_stat("PTGES", "paired_wilcox_FDR")
  ),
  expected = c(
    0.3884300540,
    0.2463247228,
    0.0299739736,
    0.0432380235,
    0.0432380235,
    0.2084128037
  ),
  tolerance = c(
    0.002,
    0.002,
    0.002,
    0.002,
    0.002,
    0.002
  )
) %>%
  dplyr::mutate(
    abs_diff = abs(observed - expected),
    pass = abs_diff <= tolerance
  )

readr::write_csv(
  qc_table,
  file.path(outdir, "FigS23_QC_check.csv")
)

cat("\n==================== S23 QC CHECK ====================\n")
print(qc_table)

if (all(qc_table$pass)) {
  cat("\nS23 QC RESULT: PASS. Values match the locked final caption/main-text logic.\n")
} else {
  cat("\nS23 QC RESULT: FAIL. Some values do not match the locked final caption/main-text logic.\n")
  stop("S23 QC failed. Please inspect FigS23_QC_check.csv and FigS23_patient_level_stats.csv.")
}

# -------------------------
# Final plot
# No in-figure 'Supplementary Figure S23' title
# -------------------------

s23_effects <- s23_effects %>%
  dplyr::mutate(
    gene = factor(gene, levels = c("SPP1", "MIF", "PTGES")),
    patient = factor(as.character(patient), levels = sort(unique(as.character(patient))))
  )

label_df <- s23_stats %>%
  dplyr::mutate(
    gene = factor(gene, levels = c("SPP1", "MIF", "PTGES")),
    label = paste0("FDR = ", signif(paired_wilcox_FDR, 3))
  )

annot_df <- s23_effects %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    x = 1,
    y = max(effect_High_minus_Low, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::left_join(label_df, by = "gene")

p_s23 <- ggplot2::ggplot(
  s23_effects,
  ggplot2::aes(x = patient, y = effect_High_minus_Low)
) +
  ggplot2::geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey40"
  ) +
  ggplot2::geom_col(
    fill = "#5DA5DA",
    color = "black",
    width = 0.75
  ) +
  ggplot2::facet_wrap(~ gene, scales = "free_y", nrow = 1) +
  ggplot2::geom_text(
    data = annot_df,
    ggplot2::aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = -0.3,
    size = 3.5
  ) +
  ggplot2::labs(
    title = NULL,
    subtitle = NULL,
    x = "Patient",
    y = "Mean expression effect (High - Low)"
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(face = "bold"),
    axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1)
  )

ggplot2::ggsave(
  file.path(outdir, "FigS23_patient_level_ligand_effects.pdf"),
  p_s23,
  width = 10.0,
  height = 4.2
)

ggplot2::ggsave(
  file.path(outdir, "FigS23_patient_level_ligand_effects.png"),
  p_s23,
  width = 10.0,
  height = 4.2,
  dpi = 300
)

cat("\nSaved S23 files:\n")
cat(file.path(outdir, "FigS23_patient_level_ligand_effects.pdf"), "\n")
cat(file.path(outdir, "FigS23_patient_level_ligand_effects.png"), "\n")
cat(file.path(outdir, "FigS23_patient_level_means.csv"), "\n")
cat(file.path(outdir, "FigS23_patient_level_effects.csv"), "\n")
cat(file.path(outdir, "FigS23_patient_level_stats.csv"), "\n")
cat(file.path(outdir, "FigS23_QC_check.csv"), "\n")
