# Supplementary Figure S24: TCGA-LIHC and GSE14520 3-year OS calibration
# Manual data-frame version copied from existing calibration summary text files.
# Output:
#   FigS24_calibration_TCGA_GSE14520_vector_clean.pdf
#   FigS24_calibration_check.txt

source("00_config_paths.R")
ensure_packages(c("ggplot2", "cowplot"))
suppressPackageStartupMessages({library(ggplot2); library(cowplot)})

tcga_cal <- data.frame(cohort = "TCGA-LIHC", row_id = 1:5,
  mean.predicted = c(0.3635342, 0.5921683, 0.6965499, 0.7579478, 0.8222000),
  KM = c(0.3268410, 0.5448811, 0.6740546, 0.8580009, 0.8183085),
  KM.corrected = c(0.3429995, 0.5502722, 0.6793886, 0.8517799, 0.8106478),
  std.err = c(0.21208150, 0.13830117, 0.10091370, 0.05550317, 0.06955676),
  Lower = c(-0.17804833, -0.20910428, -0.14725268, -0.02191906, -0.10681511),
  Upper = c(0.1347479, 0.1303588, 0.1220971, 0.2275887, 0.1376938))

gse_cal <- data.frame(cohort = "GSE14520", row_id = 1:4,
  mean.predicted = c(0.4664987, 0.6487969, 0.7209262, 0.8367237),
  KM = c(0.4361771, 0.5950716, 0.8336192, 0.8095238),
  KM.corrected = c(0.4607791, 0.5930942, 0.8125020, 0.7932815),
  std.err = c(0.16400614, 0.11817156, 0.06075358, 0.06701257),
  Lower = c(-0.1396346, -0.2172149, -0.0746497, -0.1253029),
  Upper = c(0.14492620, 0.08575479, 0.27607915, 0.07747712))

tcga_cal$ymin <- pmax(0, tcga_cal$KM + tcga_cal$Lower); tcga_cal$ymax <- pmin(1, tcga_cal$KM + tcga_cal$Upper)
gse_cal$ymin <- pmax(0, gse_cal$KM + gse_cal$Lower); gse_cal$ymax <- pmin(1, gse_cal$KM + gse_cal$Upper)

check_file <- "FigS24_calibration_check.txt"
sink(check_file)
cat("===== FigS24 calibration check =====\n\n")
cat("TCGA-LIHC 3-year OS calibration\nTime point: 1095 days\nBootstrap B = 100\nGroup size m = 60\nn = 341\n\n"); print(tcga_cal)
cat("\n\nGSE14520 3-year OS calibration\nTime point: 36 months\nBootstrap B = 100\nGroup size m = 45\nn = 217\n\n"); print(gse_cal)
sink()

make_cal_plot <- function(df, title_text, xlim_use, ylim_use) {
  plot_df <- df[order(df$mean.predicted), ]
  ggplot(plot_df, aes(x = mean.predicted)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", linewidth = 0.6) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.012, linewidth = 0.55) +
    geom_line(aes(y = KM), linewidth = 0.75) + geom_point(aes(y = KM), size = 2.2) +
    geom_line(aes(y = KM.corrected), linetype = "dotted", linewidth = 0.85) + geom_point(aes(y = KM.corrected), size = 1.8) +
    coord_cartesian(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
    labs(title = title_text, x = "Predicted 3-year overall survival probability", y = "Observed 3-year overall survival probability") +
    theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", hjust = 0.5), axis.title = element_text(size = 11), axis.text = element_text(size = 10), plot.margin = margin(8, 10, 8, 8))
}
p_final <- cowplot::plot_grid(make_cal_plot(tcga_cal, "TCGA-LIHC", c(0.30, 0.86), c(0.20, 1.00)), make_cal_plot(gse_cal, "GSE14520", c(0.42, 0.88), c(0.30, 1.00)), labels = c("A", "B"), nrow = 1, label_size = 14, rel_widths = c(1, 1))
ggsave("FigS24_calibration_TCGA_GSE14520_vector_clean.pdf", p_final, device = cairo_pdf, width = 9.2, height = 4.6, units = "in")
cat("Saved figure and check file.\n")
