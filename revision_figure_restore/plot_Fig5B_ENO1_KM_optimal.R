# Figure 5B: ENO1 Kaplan-Meier survival curve, optimal cutoff
# Output:
#   Fig5B_ENO1_KM_optimal_vector_clean.pdf
#   Fig5B_ENO1_KM_optimal_check.txt

source("00_config_paths.R")
ensure_packages(c("readr", "survival", "survminer", "ggplot2", "cowplot", "AnnotationDbi", "org.Hs.eg.db"))
suppressPackageStartupMessages({
  library(readr); library(survival); library(survminer); library(ggplot2); library(cowplot)
  library(AnnotationDbi); library(org.Hs.eg.db)
})

expr_file <- "TCGA-LIHC.star_fpkm.tsv.gz"
surv_file <- "TCGA-LIHC.survival.tsv.gz"
if (!file.exists(expr_file)) stop("Missing expression file: ", expr_file)
if (!file.exists(surv_file)) stop("Missing survival file: ", surv_file)

fpkm <- readr::read_tsv(expr_file, show_col_types = FALSE)
gene_id <- as.character(fpkm[[1]])
expr_mat <- as.matrix(fpkm[, -1]); mode(expr_mat) <- "numeric"; rownames(expr_mat) <- gene_id
ensembl_ids <- gsub("\\..*", "", rownames(expr_mat))
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
gene_symbols <- as.character(gene_symbols)
rownames(expr_mat) <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)
if (any(duplicated(rownames(expr_mat)))) {
  gene_means <- rowMeans(expr_mat, na.rm = TRUE)
  keep_idx <- unlist(tapply(seq_along(gene_means), rownames(expr_mat), function(ii) ii[which.max(gene_means[ii])]))
  expr_mat <- expr_mat[keep_idx, , drop = FALSE]
}
if (!"ENO1" %in% rownames(expr_mat)) stop("ENO1 not found in expression matrix.")
tumor_samples <- colnames(expr_mat)[substr(colnames(expr_mat), 14, 15) == "01"]
expr_df <- data.frame(sample_full = tumor_samples,
                      sample16 = substr(tumor_samples, 1, 16),
                      ENO1_expr = as.numeric(log2(expr_mat["ENO1", tumor_samples] + 1)),
                      stringsAsFactors = FALSE)
expr_df <- expr_df[!duplicated(expr_df$sample16), ]

surv_raw <- as.data.frame(readr::read_tsv(surv_file, show_col_types = FALSE))
if (!all(c("sample", "OS", "OS.time") %in% colnames(surv_raw))) stop("Survival file missing sample/OS/OS.time.")
surv_df <- data.frame(sample16 = substr(as.character(surv_raw$sample), 1, 16),
                      OS = as.numeric(surv_raw$OS), OS.time = as.numeric(surv_raw$OS.time), stringsAsFactors = FALSE)
surv_df <- surv_df[!duplicated(surv_df$sample16), ]
plot_df <- merge(expr_df, surv_df, by = "sample16")
plot_df <- plot_df[complete.cases(plot_df[, c("ENO1_expr", "OS", "OS.time")]), ]

cut_result <- surv_cutpoint(plot_df, time = "OS.time", event = "OS", variables = "ENO1_expr", minprop = 0.1)
cutoff <- cut_result$cutpoint$cutpoint[1]
plot_df$ENO1_group <- ifelse(plot_df$ENO1_expr > cutoff, "ENO1-high", "ENO1-low")
plot_df$ENO1_group <- factor(plot_df$ENO1_group, levels = c("ENO1-low", "ENO1-high"))
fit <- survfit(Surv(OS.time, OS) ~ ENO1_group, data = plot_df)
lr <- survdiff(Surv(OS.time, OS) ~ ENO1_group, data = plot_df)
pval <- 1 - pchisq(lr$chisq, df = 1)

check_file <- "Fig5B_ENO1_KM_optimal_check.txt"
sink(check_file)
cat("===== Fig5B ENO1 KM optimal cutoff check =====\n\n")
cat("Matched samples with OS:", nrow(plot_df), "\n")
cat("Events:", sum(plot_df$OS == 1), "\n")
cat("Optimal cutoff:", cutoff, "\n")
cat("Group counts:\n"); print(table(plot_df$ENO1_group))
cat("Log-rank p:", pval, "\n")
sink()

g <- ggsurvplot(fit, data = plot_df, pval = TRUE, risk.table = TRUE, risk.table.height = 0.25,
                conf.int = FALSE, censor = TRUE, palette = c("#4DBBD5", "#E64B35"),
                legend.title = "", legend.labs = c("ENO1-low", "ENO1-high"),
                xlab = "Time (days)", ylab = "Overall survival probability",
                ggtheme = theme_classic(base_size = 13), tables.theme = theme_classic(base_size = 12))
km_plot <- g$plot + theme(legend.position = "top", plot.margin = margin(5, 10, 0, 10))
risk_table <- g$table + theme(legend.position = "none", plot.margin = margin(0, 10, 5, 10))
combined <- cowplot::plot_grid(km_plot, risk_table, ncol = 1, align = "v", axis = "lr", rel_heights = c(3.2, 1.05))
ggsave("Fig5B_ENO1_KM_optimal_vector_clean.pdf", combined, device = cairo_pdf, width = 7.2, height = 6.6, units = "in")
cat("Saved figure and check file.\n")
