# ============================================================
# 08_glycolysis_gradient.R
#
# Purpose:
#   Reproduce Supplementary Figure S15:
#   immunosuppressive ligand expression along the continuous
#   glycolysis activity gradient in tumor-derived hepatocytes.
#
# Manuscript sections:
#   Methods Section 2.11
#   Results Section 3.12
#   Supplementary Figure S15
#
# Final restored method:
#   - Input object: seurat_final.rds
#   - Cells: tumor-derived hepatocytes
#   - Required cell count: 15,391
#   - Glycolysis column: Glycolysis_AUC
#   - Binning: 10 equal-width bins spanning the Glycolysis_AUC range
#   - Expression: mean log-normalized RNA expression
#   - Genes: ENO1, LDHA, SPP1, MIF, PTGES
#
# Output:
#   - FigS15_source_data.csv
#   - FigS15_QC_check.csv
#   - FigS15_final_no_title.pdf
#   - FigS15_final_no_title.png
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(scales)
})

# ---------- avoid select() masking ----------
select     <- dplyr::select
mutate     <- dplyr::mutate
filter     <- dplyr::filter
arrange    <- dplyr::arrange
summarise  <- dplyr::summarise
group_by   <- dplyr::group_by
ungroup    <- dplyr::ungroup
left_join  <- dplyr::left_join
rename     <- dplyr::rename
count      <- dplyr::count
bind_rows  <- dplyr::bind_rows
slice      <- dplyr::slice
slice_head <- dplyr::slice_head
n_distinct <- dplyr::n_distinct
all_of     <- dplyr::all_of

# ============================================================
# User-configurable paths
# ============================================================

# Recommended use:
#   Sys.setenv(SEURAT_RDS = "D:/scRNA_project/seurat_final.rds")
#   Sys.setenv(FIGS15_OUTDIR = "D:/scRNA_project/rerun_S12_S15_S23_final")
#   source("08_glycolysis_gradient.R")

candidate_rds_paths <- c(
  Sys.getenv("SEURAT_RDS", unset = NA_character_),
  "seurat_final.rds",
  "D:/scRNA_project/seurat_final.rds"
)

candidate_rds_paths <- candidate_rds_paths[!is.na(candidate_rds_paths)]
candidate_rds_paths <- candidate_rds_paths[file.exists(candidate_rds_paths)]

if (length(candidate_rds_paths) == 0) {
  stop(
    "Cannot find seurat_final.rds. Set the path first, for example:\n",
    "Sys.setenv(SEURAT_RDS = 'D:/scRNA_project/seurat_final.rds')"
  )
}

rds_path <- candidate_rds_paths[1]

outdir <- Sys.getenv(
  "FIGS15_OUTDIR",
  unset = "rerun_S12_S15_S23_final"
)

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("Using RDS:\n")
print(rds_path)

cat("Using output directory:\n")
print(normalizePath(outdir, winslash = "/", mustWork = FALSE))

# ============================================================
# Fixed metadata columns
# ============================================================

patient_col  <- "patient"
celltype_col <- "cell_type"
site_col     <- "site"
gly_col      <- "Glycolysis_AUC"

required_meta_cols <- c(patient_col, celltype_col, site_col, gly_col)

# ============================================================
# Load object
# ============================================================

seu <- readRDS(rds_path)

if (!inherits(seu, "Seurat")) {
  stop("Input object is not a Seurat object.")
}

meta <- seu@meta.data %>%
  tibble::rownames_to_column("cell")

missing_meta_cols <- setdiff(required_meta_cols, colnames(meta))

if (length(missing_meta_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_meta_cols, collapse = ", ")
  )
}

cat("\nObject class:\n")
print(class(seu))

cat("\nTotal cells:\n")
print(ncol(seu))

# ============================================================
# Normalized expression matrix
# ============================================================

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

cat("\nAssay used:\n")
print(assay_use)

cat("\nNormalized matrix dimensions:\n")
print(dim(norm_mat))

# ============================================================
# Select tumor-derived hepatocytes
# ============================================================

tumor_hep <- meta %>%
  dplyr::filter(
    .data[[site_col]] == "Tumor",
    .data[[celltype_col]] %in% c("Hepatocyte")
  ) %>%
  dplyr::mutate(
    Glycolysis_AUC = .data[[gly_col]],
    patient = .data[[patient_col]]
  )

cat("\nSelected tumor-derived hepatocytes:\n")
print(nrow(tumor_hep))

if (nrow(tumor_hep) != 15391) {
  stop(
    "Selected tumor-derived hepatocytes != 15391. ",
    "Check cell_type/site labels before plotting."
  )
}

median_cut <- median(tumor_hep$Glycolysis_AUC, na.rm = TRUE)

tumor_hep <- tumor_hep %>%
  dplyr::mutate(
    gly_group = ifelse(Glycolysis_AUC > median_cut, "High", "Low"),
    gly_group = factor(gly_group, levels = c("Low", "High"))
  )

cat("\nMedian Glycolysis_AUC cutoff:\n")
print(median_cut)

cat("\nGlyco group counts:\n")
print(table(tumor_hep$gly_group))

if (as.integer(table(tumor_hep$gly_group)["Low"]) != 7696 ||
    as.integer(table(tumor_hep$gly_group)["High"]) != 7695) {
  stop("Glyco group counts do not match expected Low=7696 / High=7695.")
}

# ============================================================
# Expression helper
# ============================================================

pull_gene_expr <- function(genes, cells, norm_mat) {
  genes_present <- intersect(genes, rownames(norm_mat))

  if (length(genes_present) != length(genes)) {
    warning(
      "Missing genes: ",
      paste(setdiff(genes, genes_present), collapse = ", ")
    )
  }

  expr <- as.matrix(norm_mat[genes_present, cells, drop = FALSE])
  expr <- t(expr)
  expr <- as.data.frame(expr)
  expr <- tibble::rownames_to_column(expr, "cell")
  expr
}

# ============================================================
# Supplementary Figure S15 source data
# ============================================================

source_csv <- file.path(outdir, "FigS15_source_data.csv")
qc_csv     <- file.path(outdir, "FigS15_QC_check.csv")
pdf_out    <- file.path(outdir, "FigS15_final_no_title.pdf")
png_out    <- file.path(outdir, "FigS15_final_no_title.png")

s15_genes <- c("ENO1", "LDHA", "SPP1", "MIF", "PTGES")

gene_class_df <- tibble::tibble(
  gene = s15_genes,
  gene_class = c(
    "Glycolytic enzyme",
    "Glycolytic enzyme",
    "Immunosuppressive ligand",
    "Immunosuppressive ligand",
    "Immunosuppressive ligand"
  ),
  line_type = c("solid", "solid", "dashed", "dashed", "dashed")
)

genes_missing <- setdiff(s15_genes, rownames(norm_mat))

if (length(genes_missing) > 0) {
  stop("Missing genes in normalized matrix: ", paste(genes_missing, collapse = ", "))
}

# Important:
# This is equal-width binning across the Glycolysis_AUC range.
# It is NOT equal-cell-count ntile() binning.
s15_cells <- tumor_hep %>%
  dplyr::select(cell, patient, Glycolysis_AUC, gly_group) %>%
  dplyr::arrange(Glycolysis_AUC, cell) %>%
  dplyr::mutate(
    gly_bin = as.integer(cut(
      Glycolysis_AUC,
      breaks = 10,
      include.lowest = TRUE,
      labels = FALSE
    )),
    gly_bin_label = paste0("Bin ", gly_bin)
  )

cat("\nS15 equal-width bin counts:\n")
print(table(s15_cells$gly_bin))

if (nrow(s15_cells) != 15391) {
  stop("S15 selected cells != 15391.")
}

if (length(unique(s15_cells$gly_bin)) != 10) {
  stop("S15 bin count != 10.")
}

s15_expr_wide <- pull_gene_expr(
  genes = s15_genes,
  cells = s15_cells$cell,
  norm_mat = norm_mat
)

s15_expr_long <- s15_expr_wide %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(s15_genes),
    names_to = "gene",
    values_to = "expr"
  ) %>%
  dplyr::left_join(
    s15_cells %>%
      dplyr::select(cell, patient, Glycolysis_AUC, gly_group, gly_bin, gly_bin_label),
    by = "cell"
  ) %>%
  dplyr::left_join(gene_class_df, by = "gene")

s15_source <- s15_expr_long %>%
  dplyr::group_by(gly_bin, gly_bin_label, gene, gene_class, line_type) %>%
  dplyr::summarise(
    n_cells = dplyr::n_distinct(cell),
    n_patients = dplyr::n_distinct(patient),
    gly_min = min(Glycolysis_AUC, na.rm = TRUE),
    gly_q25 = as.numeric(stats::quantile(Glycolysis_AUC, 0.25, na.rm = TRUE)),
    gly_mean = mean(Glycolysis_AUC, na.rm = TRUE),
    gly_median = median(Glycolysis_AUC, na.rm = TRUE),
    gly_q75 = as.numeric(stats::quantile(Glycolysis_AUC, 0.75, na.rm = TRUE)),
    gly_max = max(Glycolysis_AUC, na.rm = TRUE),
    mean_expr = mean(expr, na.rm = TRUE),
    median_expr = median(expr, na.rm = TRUE),
    pct_detected = mean(expr > 0, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  dplyr::arrange(gene, gly_bin)

readr::write_csv(s15_source, source_csv)

cat("\nS15 source CSV written:\n")
print(source_csv)

# ============================================================
# QC
# ============================================================

get_mean <- function(gene_name, bin_id) {
  s15_source %>%
    dplyr::filter(gene == gene_name, gly_bin == bin_id) %>%
    dplyr::pull(mean_expr)
}

peak_bin <- function(gene_name) {
  s15_source %>%
    dplyr::filter(gene == gene_name) %>%
    dplyr::arrange(dplyr::desc(mean_expr), gly_bin) %>%
    dplyr::slice(1) %>%
    dplyr::pull(gly_bin)
}

round3_equal <- function(x, expected) {
  isTRUE(round(as.numeric(x), 3) == expected)
}

bin_counts <- s15_cells %>%
  dplyr::count(gly_bin, name = "n_cells") %>%
  dplyr::arrange(gly_bin)

gly_mean_by_bin <- s15_cells %>%
  dplyr::group_by(gly_bin) %>%
  dplyr::summarise(gly_mean = mean(Glycolysis_AUC, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(gly_bin) %>%
  dplyr::pull(gly_mean)

spp1_bins_7_9 <- s15_source %>%
  dplyr::filter(gene == "SPP1", gly_bin %in% 7:9) %>%
  dplyr::arrange(gly_bin)

mif_values <- s15_source %>%
  dplyr::filter(gene == "MIF") %>%
  dplyr::arrange(gly_bin)

ptges_values <- s15_source %>%
  dplyr::filter(gene == "PTGES") %>%
  dplyr::arrange(gly_bin)

s15_qc <- tibble::tibble(
  check_id = c(
    "selected_tumor_hepatocytes_n",
    "glyco_low_count",
    "glyco_high_count",
    "s15_gene_count",
    "s15_missing_genes",
    "s15_bin_method",
    "s15_expression_method",
    "s15_bin_count",
    "s15_total_cells_by_bin",
    "glycolysis_bin_mean_increases",
    "ENO1_bin1_matches_caption_round3",
    "ENO1_bin10_matches_caption_round3",
    "LDHA_bin1_matches_caption_round3",
    "LDHA_bin10_matches_caption_round3",
    "SPP1_peak_bin_7_to_9",
    "SPP1_bins_7_to_9_match_caption_range",
    "MIF_bin8_gt_bin1",
    "PTGES_recorded_not_hard_block"
  ),
  observed = c(
    as.character(nrow(tumor_hep)),
    as.character(as.integer(table(tumor_hep$gly_group)["Low"])),
    as.character(as.integer(table(tumor_hep$gly_group)["High"])),
    as.character(length(unique(s15_source$gene))),
    ifelse(length(genes_missing) == 0, "none", paste(genes_missing, collapse = ";")),
    "equal-width cut(Glycolysis_AUC, breaks = 10)",
    "mean log-normalized expression",
    as.character(length(unique(s15_source$gly_bin))),
    as.character(sum(bin_counts$n_cells)),
    paste(round(gly_mean_by_bin, 6), collapse = " -> "),
    as.character(round(get_mean("ENO1", 1), 3)),
    as.character(round(get_mean("ENO1", 10), 3)),
    as.character(round(get_mean("LDHA", 1), 3)),
    as.character(round(get_mean("LDHA", 10), 3)),
    as.character(peak_bin("SPP1")),
    paste(round(spp1_bins_7_9$mean_expr, 3), collapse = ", "),
    as.character(get_mean("MIF", 8) > get_mean("MIF", 1)),
    paste(round(ptges_values$mean_expr, 3), collapse = " -> ")
  ),
  expected = c(
    "15391",
    "7696",
    "7695",
    "5",
    "none",
    "equal-width bins across AUCell glycolysis score range",
    "mean log-normalized expression",
    "10",
    "15391",
    "strictly increasing",
    "0.367",
    "3.283; display may show 3.28",
    "0.172",
    "3.451; display may show 3.45",
    "7, 8, or 9",
    "2.3 to 2.6",
    "TRUE",
    "record only; PTGES is secondary and not a hard blocking endpoint"
  ),
  pass = c(
    nrow(tumor_hep) == 15391,
    as.integer(table(tumor_hep$gly_group)["Low"]) == 7696,
    as.integer(table(tumor_hep$gly_group)["High"]) == 7695,
    length(unique(s15_source$gene)) == 5,
    length(genes_missing) == 0,
    TRUE,
    TRUE,
    length(unique(s15_source$gly_bin)) == 10,
    sum(bin_counts$n_cells) == 15391,
    all(diff(gly_mean_by_bin) > 0),
    round3_equal(get_mean("ENO1", 1), 0.367),
    abs(round(get_mean("ENO1", 10), 3) - 3.283) <= 0.005,
    round3_equal(get_mean("LDHA", 1), 0.172),
    abs(round(get_mean("LDHA", 10), 3) - 3.451) <= 0.005,
    peak_bin("SPP1") %in% 7:9,
    all(round(spp1_bins_7_9$mean_expr, 2) >= 2.29 &
          round(spp1_bins_7_9$mean_expr, 2) <= 2.61),
    get_mean("MIF", 8) > get_mean("MIF", 1),
    TRUE
  )
)

readr::write_csv(s15_qc, qc_csv)

cat("\nS15 QC CSV written:\n")
print(qc_csv)

cat("\nS15 QC table:\n")
print(s15_qc, n = Inf)

if (!all(s15_qc$pass)) {
  cat("\nS15 QC FAILED. Failed rows:\n")
  print(s15_qc %>% dplyr::filter(!pass), n = Inf)
  stop("Do not plot. QC failed.")
}

cat("\nS15 QC PASS. Writing PDF/PNG.\n")

# ============================================================
# Plot
# ============================================================

plot_df <- s15_source %>%
  dplyr::mutate(
    gene = factor(gene, levels = s15_genes),
    gene_class = factor(
      gene_class,
      levels = c("Glycolytic enzyme", "Immunosuppressive ligand")
    )
  )

p_s15 <- ggplot(
  plot_df,
  aes(
    x = gly_bin,
    y = mean_expr,
    color = gene,
    linetype = gene_class,
    group = gene
  )
) +
  geom_line(linewidth = 1.05) +
  geom_point(size = 2.2, stroke = 0.25) +
  scale_x_continuous(
    breaks = 1:10,
    labels = 1:10,
    expand = expansion(mult = c(0.015, 0.035))
  ) +
  scale_linetype_manual(
    values = c(
      "Glycolytic enzyme" = "solid",
      "Immunosuppressive ligand" = "dashed"
    )
  ) +
  labs(
    x = "Glycolysis activity bin",
    y = "Mean log-normalized expression",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_blank(),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.45, color = "black")
  )

tryCatch(
  {
    ggsave(pdf_out, p_s15, width = 6.8, height = 4.6, device = grDevices::cairo_pdf)
  },
  error = function(e) {
    message("cairo_pdf failed; falling back to standard pdf device.")
    ggsave(pdf_out, p_s15, width = 6.8, height = 4.6, device = "pdf")
  }
)

ggsave(png_out, p_s15, width = 6.8, height = 4.6, dpi = 600)

cat("\nS15 final files written:\n")
print(source_csv)
print(qc_csv)
print(pdf_out)
print(png_out)

cat("\nDONE: Supplementary Figure S15 reproduced.\n")
