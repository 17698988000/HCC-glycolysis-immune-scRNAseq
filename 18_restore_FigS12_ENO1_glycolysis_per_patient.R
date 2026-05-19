# ============================================================
# Supplementary Figure S12
# Per-patient ENO1-glycolysis correlation in tumor-derived hepatocytes
#
# Input:
#   seurat_final.rds
#
# Output:
#   revision_figure_restore/FigS12/FigS12_source_data.csv
#   revision_figure_restore/FigS12/FigS12_QC_check.csv
#   revision_figure_restore/FigS12/FigS12_final_two_color_no_stars.pdf
#   revision_figure_restore/FigS12/FigS12_final_two_color_no_stars.png
#
# Final current-object values:
#   selected tumor-derived hepatocytes = 15,391
#   GlycoLow = 7,696
#   GlycoHigh = 7,695
#   global Pearson R = 0.57
#   rho range = 0.363 (HCC04) to 0.633 (HCC10)
#   median rho = 0.498
#   raw p range = 4.43e-248 to 5.06e-17
#   all BH-adjusted p < 0.05
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(stringr)
})

# ---------- user inputs ----------
rds_path <- "seurat_final.rds"

patient_col  <- "patient"
celltype_col <- "cell_type"
site_col     <- "site"
gly_col      <- "Glycolysis_AUC"

outdir <- file.path("revision_figure_restore", "FigS12")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

out_source <- file.path(outdir, "FigS12_source_data.csv")
out_qc     <- file.path(outdir, "FigS12_QC_check.csv")
out_pdf    <- file.path(outdir, "FigS12_final_two_color_no_stars.pdf")
out_png    <- file.path(outdir, "FigS12_final_two_color_no_stars.png")

# ---------- load ----------
seu <- readRDS(rds_path)
meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell")

# ---------- choose assay / normalized layer ----------
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

# ---------- select tumor-derived hepatocytes ----------
tumor_hepatocyte_labels <- c("Hepatocyte")

tumor_hep <- meta %>%
  dplyr::filter(
    .data[[site_col]] == "Tumor",
    .data[[celltype_col]] %in% tumor_hepatocyte_labels
  ) %>%
  dplyr::mutate(
    Glycolysis_AUC = .data[[gly_col]],
    patient = .data[[patient_col]]
  )

cat("\nSelected tumor-derived hepatocytes:", nrow(tumor_hep), "\n")

if (nrow(tumor_hep) != 15391) {
  stop("Selected tumor-derived hepatocytes != 15391. Stop before plotting.")
}

# ---------- GlycoHigh / GlycoLow check ----------
median_cut <- median(tumor_hep$Glycolysis_AUC, na.rm = TRUE)

tumor_hep <- tumor_hep %>%
  dplyr::mutate(
    gly_group = ifelse(Glycolysis_AUC > median_cut, "High", "Low"),
    gly_group = factor(gly_group, levels = c("Low", "High"))
  )

gly_tab <- table(tumor_hep$gly_group)

cat("\nGlyco group counts:\n")
print(gly_tab)

if (as.integer(gly_tab["Low"]) != 7696 || as.integer(gly_tab["High"]) != 7695) {
  stop("Glyco group counts do not match expected Low=7696 / High=7695.")
}

# ---------- ENO1 expression ----------
if (!"ENO1" %in% rownames(norm_mat)) {
  stop("ENO1 not found in normalized expression matrix.")
}

s12_cell_df <- tumor_hep %>%
  dplyr::select(cell, patient, Glycolysis_AUC) %>%
  dplyr::mutate(
    ENO1_expr = as.numeric(norm_mat["ENO1", cell, drop = TRUE])
  ) %>%
  dplyr::filter(
    !is.na(patient),
    !is.na(Glycolysis_AUC),
    !is.na(ENO1_expr)
  )

if (nrow(s12_cell_df) != 15391) {
  stop("S12 dataframe after ENO1 merge != 15391.")
}

# ---------- global Pearson reference line ----------
global_pearson_R <- cor(
  s12_cell_df$ENO1_expr,
  s12_cell_df$Glycolysis_AUC,
  method = "pearson",
  use = "complete.obs"
)

cat("\nGlobal Pearson R:", global_pearson_R, "\n")

if (abs(global_pearson_R - 0.57) > 0.01) {
  stop("Global Pearson R is not close to expected 0.57. Stop before plotting.")
}

# ---------- per-patient Spearman ----------
calc_patient_spearman <- function(df) {
  ct <- suppressWarnings(
    cor.test(
      df$ENO1_expr,
      df$Glycolysis_AUC,
      method = "spearman",
      exact = FALSE
    )
  )

  rho <- unname(ct$estimate)
  n <- nrow(df)

  # Fisher z-transformation CI
  z <- atanh(rho)
  se <- 1 / sqrt(n - 3)
  zcrit <- qnorm(0.975)

  tibble::tibble(
    n_cells = n,
    rho = rho,
    p_value = ct$p.value,
    ci_low = tanh(z - zcrit * se),
    ci_high = tanh(z + zcrit * se)
  )
}

s12_source <- s12_cell_df %>%
  dplyr::group_by(patient) %>%
  dplyr::group_modify(~ calc_patient_spearman(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    fdr = p.adjust(p_value, method = "BH"),
    global_pearson_R = global_pearson_R
  ) %>%
  dplyr::arrange(rho) %>%
  dplyr::mutate(
    patient = factor(patient, levels = patient)
  )

readr::write_csv(s12_source, out_source)

# ---------- QC ----------
qc <- tibble::tibble(
  metric = c(
    "selected_tumor_derived_hepatocytes",
    "glyco_low_n",
    "glyco_high_n",
    "n_patients",
    "global_pearson_R",
    "rho_min",
    "rho_min_patient",
    "rho_max",
    "rho_max_patient",
    "rho_median",
    "raw_p_min",
    "raw_p_max",
    "all_FDR_lt_0.05"
  ),
  observed = c(
    as.character(nrow(tumor_hep)),
    as.character(as.integer(gly_tab["Low"])),
    as.character(as.integer(gly_tab["High"])),
    as.character(dplyr::n_distinct(s12_source$patient)),
    sprintf("%.6f", global_pearson_R),
    sprintf("%.6f", min(s12_source$rho)),
    as.character(s12_source$patient[which.min(s12_source$rho)]),
    sprintf("%.6f", max(s12_source$rho)),
    as.character(s12_source$patient[which.max(s12_source$rho)]),
    sprintf("%.6f", median(s12_source$rho)),
    format(min(s12_source$p_value), scientific = TRUE, digits = 6),
    format(max(s12_source$p_value), scientific = TRUE, digits = 6),
    as.character(all(s12_source$fdr < 0.05))
  ),
  expected = c(
    "15391",
    "7696",
    "7695",
    "8",
    "approximately 0.57",
    "approximately 0.363",
    "HCC04",
    "approximately 0.633",
    "HCC10",
    "approximately 0.498",
    "approximately 4.43e-248",
    "approximately 5.06e-17",
    "TRUE"
  )
)

readr::write_csv(qc, out_qc)

cat("\nS12 QC summary:\n")
print(qc, n = Inf)

# ---------- final two-color no-stars plot ----------
point_col <- "#F05A3A"
err_col   <- "#39B5D8"

p <- ggplot(s12_source, aes(x = patient, y = rho, size = n_cells)) +
  geom_hline(
    yintercept = global_pearson_R,
    linetype = "dashed",
    linewidth = 0.8,
    color = "grey45"
  ) +
  geom_errorbar(
    aes(ymin = ci_low, ymax = ci_high),
    width = 0.10,
    linewidth = 1.0,
    color = err_col
  ) +
  geom_point(
    shape = 16,
    color = point_col,
    alpha = 0.95
  ) +
  annotate(
    "text",
    x = 1.0,
    y = global_pearson_R + 0.015,
    label = paste0("Global Pearson R = ", sprintf("%.2f", global_pearson_R)),
    hjust = 0,
    vjust = 0,
    size = 6.2
  ) +
  scale_size_area(max_size = 14, guide = "none") +
  scale_y_continuous(
    limits = c(0.12, 0.83),
    breaks = c(0.2, 0.4, 0.6, 0.8)
  ) +
  labs(
    title = "Per-patient ENO1-Glycolysis Correlation",
    subtitle = "Tumor-derived hepatocytes, n = 15,391",
    x = "Patient",
    y = "Per-patient Spearman rho"
  ) +
  theme_classic(base_size = 20) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 24),
    plot.subtitle = element_text(hjust = 0.5, size = 20),
    axis.title.x = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    axis.text.x = element_text(size = 18, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 18, color = "black"),
    axis.line = element_line(linewidth = 1.1, color = "black"),
    axis.ticks = element_line(linewidth = 1.1, color = "black"),
    legend.position = "none"
  )

ggsave(
  filename = out_pdf,
  plot = p,
  width = 11,
  height = 7,
  device = cairo_pdf
)

ggsave(
  filename = out_png,
  plot = p,
  width = 11,
  height = 7,
  dpi = 600
)

cat("\nS12 final files written:\n")
cat("Source:", out_source, "\n")
cat("QC    :", out_qc, "\n")
cat("PDF   :", out_pdf, "\n")
cat("PNG   :", out_png, "\n")
