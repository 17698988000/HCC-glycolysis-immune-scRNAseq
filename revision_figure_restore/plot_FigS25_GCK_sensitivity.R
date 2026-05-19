# Supplementary Figure S25: GCK sensitivity analysis
# Output:
#   FigS25_GCK_sensitivity_vector_clean.pdf
#   FigS25_GCK_sensitivity_check.txt

source("00_config_paths.R")
ensure_packages(c("ggplot2", "cowplot"))
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})
summary_file <- "JTM_GCK_sensitivity_outputs/GCK_sensitivity_summary.csv"
cell_file <- "JTM_GCK_sensitivity_outputs/GCK_sensitivity_cell_level_scores.csv"
if (!file.exists(summary_file)) stop("Missing: ", summary_file)
if (!file.exists(cell_file)) stop("Missing: ", cell_file)
sum_df <- read.csv(summary_file, check.names = FALSE)
cell_df <- read.csv(cell_file, check.names = FALSE)

n_cells <- nrow(cell_df)
detect_rate <- mean(cell_df$GCK_expression > 0, na.rm = TRUE)
mean_all <- mean(cell_df$GCK_expression, na.rm = TRUE)
median_all <- median(cell_df$GCK_expression, na.rm = TRUE)
mean_pos <- mean(cell_df$GCK_expression[cell_df$GCK_expression > 0], na.rm = TRUE)
pearson_r <- cor(cell_df$glycolysis_with_GCK, cell_df$glycolysis_no_GCK, method = "pearson")
spearman_rho <- cor(cell_df$glycolysis_with_GCK, cell_df$glycolysis_no_GCK, method = "spearman")
switched_flag <- cell_df$group_with_GCK != cell_df$group_no_GCK
n_switched <- sum(switched_flag, na.rm = TRUE)
group_concordance <- mean(!switched_flag, na.rm = TRUE)
cell_df$switched <- factor(ifelse(switched_flag, "Switched", "Unchanged"), levels = c("Unchanged", "Switched"))

check_file <- "FigS25_GCK_sensitivity_check.txt"
sink(check_file)
cat("===== FigS25 GCK sensitivity check =====\n\n")
cat("summary_file:", summary_file, "\ncell_file:", cell_file, "\n\n")
cat("n_cells =", n_cells, "\n")
cat("detect_rate =", detect_rate, "\n")
cat("detect_rate_percent =", round(detect_rate * 100, 4), "%\n")
cat("mean_all_cells =", mean_all, "\nmedian_all_cells =", median_all, "\nmean_positive_cells =", mean_pos, "\n\n")
cat("Pearson_r =", pearson_r, "\nSpearman_rho =", spearman_rho, "\n\n")
cat("group_concordance =", group_concordance, "\nn_switched_cells =", n_switched, "\n\n")
cat("Confusion table:\n"); print(table(with_GCK = cell_df$group_with_GCK, no_GCK = cell_df$group_no_GCK))
cat("\nSummary CSV contents:\n"); print(sum_df)
sink()

label_A <- paste0("Tumor hepatocytes = ", n_cells, "\nDetected in ", sprintf("%.3f", detect_rate * 100), "% of cells", "\nMedian = ", round(median_all, 3), "\nMean positive = ", sprintf("%.3f", mean_pos))
pA <- ggplot(cell_df, aes(x = GCK_expression)) +
  geom_histogram(bins = 50, fill = "#4DBBD5", color = "black", linewidth = 0.25) +
  labs(title = "GCK expression distribution", x = "GCK expression", y = "Cell count") +
  annotate("label", x = Inf, y = Inf, label = label_A, hjust = 1.03, vjust = 1.05, size = 3.6, fill = "white", label.size = 0.25) +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.margin = margin(8, 10, 8, 8)) + coord_cartesian(clip = "off")
label_B <- paste0("Pearson r = ", sprintf("%.5f", pearson_r), "\nSpearman rho = ", sprintf("%.5f", spearman_rho), "\nConcordance = ", sprintf("%.4f", group_concordance * 100), "%", "\nSwitched cells = ", n_switched, "/", n_cells)
pB <- ggplot(cell_df, aes(x = glycolysis_no_GCK, y = glycolysis_with_GCK)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_point(aes(color = switched), size = 1.0, alpha = 0.6) +
  scale_color_manual(values = c("Unchanged" = "grey40", "Switched" = "#E64B35")) +
  labs(title = "Glycolysis score with vs without GCK", x = "Glycolysis score without GCK", y = "Glycolysis score with GCK", color = "") +
  annotate("label", x = Inf, y = Inf, label = label_B, hjust = 1.03, vjust = 1.05, size = 3.6, fill = "white", label.size = 0.25) +
  theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5), legend.position = c(0.82, 0.18), legend.background = element_rect(fill = "white", color = "grey70"), plot.margin = margin(8, 10, 8, 8)) + coord_cartesian(clip = "off")
p_final <- cowplot::plot_grid(pA, pB, labels = c("A", "B"), nrow = 1, label_size = 14)
ggsave("FigS25_GCK_sensitivity_vector_clean.pdf", p_final, device = cairo_pdf, width = 10.5, height = 4.8, units = "in")
cat("Saved figure and check file.\n")
