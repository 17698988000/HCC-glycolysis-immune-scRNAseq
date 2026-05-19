# Figure 6A: Four-gene score Kaplan-Meier curve
# Output:
#   Fig6A_four_gene_score_KM_vector_clean.pdf
#   Fig6A_four_gene_score_KM_check.txt

source("00_config_paths.R")
ensure_packages(c("readr", "survival", "survminer", "ggplot2", "cowplot", "AnnotationDbi", "org.Hs.eg.db"))
suppressPackageStartupMessages({library(readr); library(survival); library(survminer); library(ggplot2); library(cowplot); library(AnnotationDbi); library(org.Hs.eg.db)})
expr_file <- "TCGA-LIHC.star_fpkm.tsv.gz"; surv_file <- "TCGA-LIHC.survival.tsv.gz"
model_genes <- c("TPI1", "ENO1", "LDHA", "SLC2A1")
coef <- c(TPI1 = 0.3041908, ENO1 = 0.9639654, LDHA = 1.3404374, SLC2A1 = 0.2424239)

fpkm <- readr::read_tsv(expr_file, show_col_types = FALSE)
expr_mat <- as.matrix(fpkm[, -1]); mode(expr_mat) <- "numeric"; rownames(expr_mat) <- as.character(fpkm[[1]])
ensembl_ids <- gsub("\\..*", "", rownames(expr_mat))
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr_mat) <- ifelse(is.na(gene_symbols), ensembl_ids, as.character(gene_symbols))
if (any(duplicated(rownames(expr_mat)))) {
  gene_means <- rowMeans(expr_mat, na.rm = TRUE)
  keep_idx <- unlist(tapply(seq_along(gene_means), rownames(expr_mat), function(ii) ii[which.max(gene_means[ii])]))
  expr_mat <- expr_mat[keep_idx, , drop = FALSE]
}
missing <- setdiff(model_genes, rownames(expr_mat)); if (length(missing) > 0) stop("Missing model genes: ", paste(missing, collapse = ", "))
tumor_samples <- colnames(expr_mat)[substr(colnames(expr_mat), 14, 15) == "01"]
expr_log2 <- log2(expr_mat[model_genes, tumor_samples, drop = FALSE] + 1)
score_vec <- colSums(t(t(expr_log2) * coef[model_genes]))
score_df <- data.frame(sample_full = tumor_samples, sample16 = substr(tumor_samples, 1, 16), four_gene_score = as.numeric(score_vec), stringsAsFactors = FALSE)
score_df <- score_df[!duplicated(score_df$sample16), ]
surv_raw <- as.data.frame(readr::read_tsv(surv_file, show_col_types = FALSE))
surv_df <- data.frame(sample16 = substr(as.character(surv_raw$sample), 1, 16), OS = as.numeric(surv_raw$OS), OS.time = as.numeric(surv_raw$OS.time), stringsAsFactors = FALSE)
surv_df <- surv_df[!duplicated(surv_df$sample16), ]
km_df <- merge(score_df, surv_df, by = "sample16")
km_df <- km_df[complete.cases(km_df[, c("four_gene_score", "OS", "OS.time")]), ]
median_cut <- median(km_df$four_gene_score, na.rm = TRUE)
km_df$score_group <- factor(ifelse(km_df$four_gene_score >= median_cut, "High score", "Low score"), levels = c("Low score", "High score"))
fit <- survfit(Surv(OS.time, OS) ~ score_group, data = km_df)
lr <- survdiff(Surv(OS.time, OS) ~ score_group, data = km_df); pval <- 1 - pchisq(lr$chisq, df = 1)

check_file <- "Fig6A_four_gene_score_KM_check.txt"
sink(check_file)
cat("===== Fig6A four-gene score KM check =====\n\n")
cat("Matched samples with OS:", nrow(km_df), "\n")
cat("Events:", sum(km_df$OS == 1), "\n")
cat("Median cutoff:", median_cut, "\n")
cat("Group counts:\n"); print(table(km_df$score_group))
cat("Log-rank p:", pval, "\n")
sink()

g <- ggsurvplot(fit, data = km_df, pval = TRUE, risk.table = TRUE, risk.table.height = 0.25, conf.int = FALSE, censor = TRUE, palette = c("#4DBBD5", "#E64B35"), legend.title = "", legend.labs = c("Low score", "High score"), xlab = "Time (days)", ylab = "Overall survival probability", ggtheme = theme_classic(base_size = 13), tables.theme = theme_classic(base_size = 12))
combined <- cowplot::plot_grid(g$plot + theme(legend.position = "top", plot.margin = margin(5, 10, 0, 10)), g$table + theme(legend.position = "none", plot.margin = margin(0, 10, 5, 10)), ncol = 1, align = "v", axis = "lr", rel_heights = c(3.2, 1.05))
ggsave("Fig6A_four_gene_score_KM_vector_clean.pdf", combined, device = cairo_pdf, width = 7.2, height = 6.6, units = "in")
cat("Saved figure and check file.\n")
