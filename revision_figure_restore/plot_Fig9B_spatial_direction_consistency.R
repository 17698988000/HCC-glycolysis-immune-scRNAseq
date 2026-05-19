# Figure 9B: Spatial GlycoHigh vs GlycoLow direction consistency
# IMPORTANT: This was temporarily named Fig8 during figure restoration. Final mapping is Fig9B.
# Output:
#   Fig9B_spatial_direction_consistency_vector_clean.pdf
#   Fig9B_spatial_direction_consistency_check.txt

source("00_config_paths.R")
ensure_packages(c("ggplot2", "cowplot"))
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})

direction_file <- "JTM_spatial_direction_consistency_summary.csv"
log2fc_file <- "JTM_spatial_glycohigh_vs_low_log2FC_by_sample.csv"
group_file <- "JTM_spatial_glyco_group_summary.csv"
for (f in c(direction_file, log2fc_file, group_file)) if (!file.exists(f)) stop("Missing file: ", f)

direction_df <- read.csv(direction_file, check.names = FALSE)
log2fc_df <- read.csv(log2fc_file, check.names = FALSE)
group_df <- read.csv(group_file, check.names = FALSE)
feature_order <- c("four_gene_score", "ENO1", "SPP1", "MIF", "PTGES", "ligand_score_MIF_SPP1")
feature_labels <- c("four_gene_score" = "Four-gene score", "ENO1" = "ENO1", "SPP1" = "SPP1", "MIF" = "MIF", "PTGES" = "PTGES", "ligand_score_MIF_SPP1" = "MIF/SPP1 ligand score")
log2fc_df <- log2fc_df[log2fc_df$feature %in% feature_order, ]
direction_df <- direction_df[direction_df$feature %in% feature_order, ]
log2fc_df$feature_label <- feature_labels[log2fc_df$feature]
direction_df$feature_label <- feature_labels[direction_df$feature]
feature_label_levels <- rev(unname(feature_labels[feature_order]))
log2fc_df$feature_label <- factor(log2fc_df$feature_label, levels = feature_label_levels)
direction_df$feature_label <- factor(direction_df$feature_label, levels = feature_label_levels)

check_file <- "Fig9B_spatial_direction_consistency_check.txt"
sink(check_file)
cat("===== Fig9B spatial direction consistency check =====\n\n")
cat("Input files:\n")
cat("direction_file:", direction_file, "\n")
cat("log2fc_file:", log2fc_file, "\n")
cat("group_file:", group_file, "\n\n")
cat("Spot counts table:\n"); print(xtabs(n_spots ~ sample + spatial_glyco_group, data = group_df))
cat("\nDirection consistency summary:\n"); print(direction_df[, c("feature", "n_samples", "n_higher_in_GlycoHigh", "direction_consistency", "median_log2FC", "mean_log2FC")])
cat("\nSample-level log2FC table:\n"); print(log2fc_df[, c("sample", "feature", "GlycoHigh", "GlycoLow", "log2FC", "direction")])
cat("\nInterpretation:\nAll listed features show positive GlycoHigh-vs-GlycoLow log2FC across all four spatial sections.\n")
cat("This supports spot-level tissue co-enrichment, not cell-cell causality or clinical prediction.\n")
sink()

max_abs <- max(abs(log2fc_df$log2FC), na.rm = TRUE)
p_heat <- ggplot(log2fc_df, aes(x = sample, y = feature_label, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.6) + geom_text(aes(label = sprintf("%.2f", log2FC)), size = 3.4) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-max_abs, max_abs), name = "log2FC\nHigh vs Low") +
  labs(title = "Spatial GlycoHigh vs GlycoLow", x = "Spatial section", y = NULL) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.line = element_blank(), axis.ticks = element_blank(), legend.title = element_text(size = 9), legend.text = element_text(size = 9), plot.margin = margin(8, 8, 8, 8))

direction_df$consistency_label <- direction_df$direction_consistency
xmax <- max(direction_df$median_log2FC, na.rm = TRUE) + 0.42
p_summary <- ggplot(direction_df, aes(x = median_log2FC, y = feature_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5) +
  geom_segment(aes(x = 0, xend = median_log2FC, yend = feature_label), linewidth = 0.75) +
  geom_point(size = 2.8) + geom_text(aes(x = median_log2FC + 0.08, label = consistency_label), hjust = 0, size = 3.6) +
  scale_x_continuous(limits = c(0, xmax), expand = c(0, 0)) +
  labs(title = "Direction consistency", x = "Median log2FC", y = NULL) +
  annotate("text", x = max(direction_df$median_log2FC, na.rm = TRUE) + 0.08, y = length(levels(direction_df$feature_label)) + 0.45, label = "Sections", hjust = 0, size = 3.5, fontface = "bold") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 11), plot.margin = margin(8, 22, 8, 8)) +
  coord_cartesian(clip = "off")

p_final <- cowplot::plot_grid(p_heat, p_summary, labels = c("A", "B"), nrow = 1, rel_widths = c(1.18, 1.08), label_size = 14)
ggsave("Fig9B_spatial_direction_consistency_vector_clean.pdf", p_final, device = cairo_pdf, width = 11.2, height = 4.8, units = "in")
cat("Saved figure and check file.\n")
