#!/usr/bin/env Rscript

# ============================================================
# 08_glycolysis_gradient.R
#
# Purpose:
#   Reproduce Supplementary Figure S15 from the locked source table.
#
# Scope:
#   15,391 tumor-derived hepatocytes divided into ten approximately
#   equal-cell-count AUCell glycolysis bins.
#
# Interpretation:
#   Descriptive gradient visualization only. This is not trajectory
#   inference, functional validation, or treatment-response evidence.
# ============================================================

required_packages <- c("readr", "dplyr", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) stop("Missing package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
})

source_csv <- file.path(
  "locked_source_data", "supplementary_figures",
  "FigS15_glycolysis_gradient_LOCKED_source.csv"
)
outdir <- file.path("results", "FigS15_glycolysis_gradient")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_csv)) stop("Locked S15 source table not found: ", source_csv, call. = FALSE)
x <- readr::read_csv(source_csv, show_col_types = FALSE)
required_cols <- c("glycolysis_bin", "n_cells", "glycolysis_mean", "gene", "mean_expression", "category", "bin_label")
missing_cols <- setdiff(required_cols, names(x))
if (length(missing_cols) > 0) stop("Locked S15 source table missing: ", paste(missing_cols, collapse = ", "), call. = FALSE)

bin_sizes <- x %>% dplyr::distinct(glycolysis_bin, n_cells) %>% dplyr::arrange(glycolysis_bin)
if (sum(bin_sizes$n_cells) != 15391L) stop("Locked S15 cell-count total is not 15,391.", call. = FALSE)
if (nrow(bin_sizes) != 10L) stop("Locked S15 must contain ten bins.", call. = FALSE)

x <- x %>%
  dplyr::mutate(
    gene = factor(gene, levels = c("ENO1", "LDHA", "SPP1", "MIF", "PTGES")),
    category = factor(category, levels = c("Glycolytic enzyme", "Immunosuppressive ligand"))
  )

p <- ggplot2::ggplot(
  x,
  ggplot2::aes(x = glycolysis_bin, y = mean_expression, group = gene, color = gene, linetype = category)
) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::scale_x_continuous(breaks = 1:10, labels = paste0("Bin ", 1:10)) +
  ggplot2::labs(
    x = "Glycolysis activity bin (low to high)",
    y = "Mean log-normalized expression",
    color = "Gene",
    linetype = NULL,
    subtitle = "Tumor-derived hepatocytes, n = 15,391; equal-cell-count AUCell bins"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

ggplot2::ggsave(file.path(outdir, "FigS15_glycolysis_gradient_LOCKED.pdf"), p, width = 8.5, height = 5.8, device = grDevices::cairo_pdf)
ggplot2::ggsave(file.path(outdir, "FigS15_glycolysis_gradient_LOCKED.png"), p, width = 8.5, height = 5.8, dpi = 300, bg = "white")
readr::write_csv(x, file.path(outdir, "FigS15_glycolysis_gradient_LOCKED_source.csv"))
readr::write_csv(bin_sizes, file.path(outdir, "FigS15_glycolysis_gradient_LOCKED_bin_sizes.csv"))

qc <- tibble::tibble(
  check = c("total_cells", "n_bins", "bin_1_cells", "bins_2_to_10_cells", "interpretation"),
  observed = c(sum(bin_sizes$n_cells), nrow(bin_sizes), bin_sizes$n_cells[1], paste(bin_sizes$n_cells[2:10], collapse = ";"), "descriptive gradient visualization"),
  expected = c(15391, 10, 1540, paste(rep(1539, 9), collapse = ";"), "descriptive gradient visualization")
)
readr::write_csv(qc, file.path(outdir, "FigS15_glycolysis_gradient_LOCKED_QC_check.csv"))
message("Supplementary Figure S15 locked redraw complete: ", outdir)
