# Supplementary Figure S20: RCTD-estimated spatial cell-type composition
# Output:
#   FigS20_RCTD_celltype_composition_vector_clean.pdf
#   FigS20_RCTD_celltype_composition_check.txt

source("00_config_paths.R")
ensure_packages(c("ggplot2", "tidyr"))
suppressPackageStartupMessages({library(ggplot2); library(tidyr)})

files <- c(
  "spatial_composition_regression/HCC1R_RCTD_celltype_proportions.csv",
  "spatial_composition_regression/HCC2R_RCTD_celltype_proportions.csv",
  "spatial_composition_regression/HCC3R_RCTD_celltype_proportions.csv",
  "spatial_composition_regression/HCC4R_RCTD_celltype_proportions.csv"
)
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) stop("Missing RCTD composition CSV: ", paste(missing_files, collapse = ", "))
rctd_df <- do.call(rbind, lapply(files, function(f) read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)))
celltypes <- c("B", "Endothelial", "Fibroblast", "Hepatocyte", "Myeloid", "T_NK")
if (!all(celltypes %in% colnames(rctd_df))) stop("Missing celltype columns: ", paste(setdiff(celltypes, colnames(rctd_df)), collapse = ", "))
rctd_sub <- rctd_df[, c("sample", "barcode", celltypes)]
rctd_long <- tidyr::pivot_longer(rctd_sub, cols = all_of(celltypes), names_to = "celltype", values_to = "proportion")
rctd_long$sample <- factor(rctd_long$sample, levels = c("HCC1R", "HCC2R", "HCC3R", "HCC4R"))
rctd_long$celltype <- factor(rctd_long$celltype, levels = c("Hepatocyte", "Myeloid", "T_NK", "Fibroblast", "Endothelial", "B"))
mean_df <- aggregate(proportion ~ sample + celltype, data = rctd_long, FUN = function(x) mean(x, na.rm = TRUE)); colnames(mean_df)[3] <- "mean_proportion"
median_df <- aggregate(proportion ~ sample + celltype, data = rctd_long, FUN = function(x) median(x, na.rm = TRUE)); colnames(median_df)[3] <- "median_proportion"
n_df <- aggregate(proportion ~ sample + celltype, data = rctd_long, FUN = length); colnames(n_df)[3] <- "n_spots"
summary_df <- merge(merge(mean_df, median_df, by = c("sample", "celltype")), n_df, by = c("sample", "celltype"))
summary_df$sample <- factor(summary_df$sample, levels = c("HCC1R", "HCC2R", "HCC3R", "HCC4R"))
summary_df$celltype <- factor(summary_df$celltype, levels = c("Hepatocyte", "Myeloid", "T_NK", "Fibroblast", "Endothelial", "B"))
summary_df <- summary_df[order(summary_df$sample, summary_df$celltype), ]

check_file <- "FigS20_RCTD_celltype_composition_check.txt"
sink(check_file)
cat("===== FigS20 RCTD cell-type composition check =====\n\n")
cat("Input files:\n"); print(files)
cat("\nMerged RCTD table dimensions:\n"); print(dim(rctd_df))
cat("\nSpot counts by sample:\n"); print(table(rctd_df$sample))
cat("\nCell-type columns:\n"); print(celltypes)
cat("\nMean and median RCTD proportions by sample and cell type:\n"); print(summary_df)
cat("\nCheck row sums of cell-type proportions:\n"); row_sum <- rowSums(rctd_df[, celltypes], na.rm = TRUE); print(summary(row_sum)); cat("\nFirst 10 row sums:\n"); print(head(row_sum, 10))
sink()

p <- ggplot(summary_df, aes(x = sample, y = mean_proportion, fill = celltype)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = function(x) paste0(round(x * 100), "%")) +
  scale_fill_manual(values = c("Hepatocyte" = "#E64B35", "Myeloid" = "#4DBBD5", "T_NK" = "#00A087", "Fibroblast" = "#3C5488", "Endothelial" = "#F39B7F", "B" = "#8491B4")) +
  labs(x = "Spatial section", y = "Mean RCTD-estimated proportion", fill = "Cell type") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size = 11), legend.text = element_text(size = 10), plot.margin = margin(10, 12, 10, 10))
ggsave("FigS20_RCTD_celltype_composition_vector_clean.pdf", p, device = cairo_pdf, width = 6.8, height = 4.8, units = "in")
cat("Saved figure and check file.\n")
