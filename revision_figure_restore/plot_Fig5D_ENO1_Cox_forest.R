# Figure 5D: ENO1 multivariate Cox forest plot
# Output:
#   Fig5D_ENO1_Cox_forest_vector_clean.pdf
#   Fig5D_ENO1_Cox_forest_check.txt

source("00_config_paths.R")
ensure_packages(c("readr", "survival", "survminer", "ggplot2", "AnnotationDbi", "org.Hs.eg.db"))
suppressPackageStartupMessages({library(readr); library(survival); library(survminer); library(ggplot2); library(AnnotationDbi); library(org.Hs.eg.db)})
expr_file <- "TCGA-LIHC.star_fpkm.tsv.gz"; surv_file <- "TCGA-LIHC.survival.tsv.gz"; clin_file <- "TCGA-LIHC.clinical.tsv.gz"

fpkm <- readr::read_tsv(expr_file, show_col_types = FALSE)
gene_id <- as.character(fpkm[[1]])
expr_mat <- as.matrix(fpkm[, -1]); mode(expr_mat) <- "numeric"; rownames(expr_mat) <- gene_id
ensembl_ids <- gsub("\\..*", "", rownames(expr_mat))
gene_symbols <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
rownames(expr_mat) <- ifelse(is.na(gene_symbols), ensembl_ids, as.character(gene_symbols))
if (any(duplicated(rownames(expr_mat)))) {
  gene_means <- rowMeans(expr_mat, na.rm = TRUE)
  keep_idx <- unlist(tapply(seq_along(gene_means), rownames(expr_mat), function(ii) ii[which.max(gene_means[ii])]))
  expr_mat <- expr_mat[keep_idx, , drop = FALSE]
}
tumor_samples <- colnames(expr_mat)[substr(colnames(expr_mat), 14, 15) == "01"]
expr_df <- data.frame(sample16 = substr(tumor_samples, 1, 16), ENO1_expr = as.numeric(log2(expr_mat["ENO1", tumor_samples] + 1)), stringsAsFactors = FALSE)
expr_df <- expr_df[!duplicated(expr_df$sample16), ]
surv_raw <- as.data.frame(readr::read_tsv(surv_file, show_col_types = FALSE))
surv_df <- data.frame(sample16 = substr(as.character(surv_raw$sample), 1, 16), OS = as.numeric(surv_raw$OS), OS.time = as.numeric(surv_raw$OS.time), stringsAsFactors = FALSE)
surv_df <- surv_df[!duplicated(surv_df$sample16), ]
clin_raw <- as.data.frame(readr::read_tsv(clin_file, show_col_types = FALSE))
stage_candidates <- colnames(clin_raw)[grepl("ajcc.*pathologic.*stage|pathologic.*stage|tumor.*stage|stage", colnames(clin_raw), ignore.case = TRUE)]
stage_col <- stage_candidates[grepl("ajcc.*pathologic.*stage", stage_candidates, ignore.case = TRUE)][1]
if (is.na(stage_col)) stage_col <- stage_candidates[grepl("pathologic.*stage", stage_candidates, ignore.case = TRUE)][1]
if (is.na(stage_col)) stage_col <- stage_candidates[1]
clin_df <- data.frame(sample16 = substr(as.character(clin_raw$sample), 1, 16),
                      age_at_index = suppressWarnings(as.numeric(clin_raw[["age_at_index.demographic"]])),
                      gender = as.character(clin_raw[["gender.demographic"]]),
                      ajcc_pathologic_stage = as.character(clin_raw[[stage_col]]), stringsAsFactors = FALSE)
clin_df <- clin_df[!duplicated(clin_df$sample16), ]
plot_df <- merge(expr_df, surv_df, by = "sample16"); plot_df <- merge(plot_df, clin_df, by = "sample16")
plot_df <- plot_df[complete.cases(plot_df[, c("ENO1_expr", "OS", "OS.time", "age_at_index", "gender", "ajcc_pathologic_stage")]), ]
cut_result <- surv_cutpoint(plot_df, time = "OS.time", event = "OS", variables = "ENO1_expr", minprop = 0.1)
cutoff <- cut_result$cutpoint$cutpoint[1]
plot_df$ENO1_group <- factor(ifelse(plot_df$ENO1_expr > cutoff, "High", "Low"), levels = c("Low", "High"))
stage_upper <- toupper(plot_df$ajcc_pathologic_stage)
plot_df$stage_binary <- ifelse(grepl("III|IV", stage_upper), "III-IV", ifelse(grepl("I|II", stage_upper), "I-II", NA))
plot_df <- plot_df[!is.na(plot_df$stage_binary), ]
plot_df$stage_binary <- factor(plot_df$stage_binary, levels = c("I-II", "III-IV"))
plot_df$gender <- factor(ifelse(toupper(plot_df$gender) %in% c("MALE", "M"), "male", "female"), levels = c("female", "male"))
cox_fit <- coxph(Surv(OS.time, OS) ~ ENO1_group + age_at_index + gender + stage_binary, data = plot_df)
cox_sum <- summary(cox_fit)
forest_df <- data.frame(variable = rownames(cox_sum$coefficients), HR = cox_sum$conf.int[, "exp(coef)"], lower = cox_sum$conf.int[, "lower .95"], upper = cox_sum$conf.int[, "upper .95"], p = cox_sum$coefficients[, "Pr(>|z|)"], stringsAsFactors = FALSE)
forest_df$label <- c("ENO1 high vs low", "Age, per year", "Male vs female", "AJCC stage III-IV vs I-II")
forest_df$p_label <- ifelse(forest_df$p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", forest_df$p)))
forest_df$hr_label <- sprintf("%.2f (%.2f-%.2f)", forest_df$HR, forest_df$lower, forest_df$upper)
forest_df$label <- factor(forest_df$label, levels = rev(forest_df$label))

check_file <- "Fig5D_ENO1_Cox_forest_check.txt"
sink(check_file)
cat("===== Fig5D ENO1 Cox forest check =====\n\n")
cat("Using stage column:", stage_col, "\n")
cat("Final Cox samples:", nrow(plot_df), "\n")
cat("Events:", sum(plot_df$OS == 1), "\n")
cat("Optimal cutoff:", cutoff, "\n")
cat("Coefficients:\n"); print(cox_sum$coefficients)
cat("Hazard ratios:\n"); print(cox_sum$conf.int)
cat("C-index:\n"); print(cox_sum$concordance)
sink()

p <- ggplot(forest_df, aes(x = HR, y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.18, linewidth = 0.7) +
  geom_point(size = 2.7) +
  geom_text(aes(x = 3.75, label = hr_label), hjust = 0, size = 3.7) +
  geom_text(aes(x = 5.4, label = p_label), hjust = 0, size = 3.7) +
  scale_x_continuous(limits = c(0.4, 6.2), breaks = c(0.5, 1, 2, 3, 4, 5, 6)) +
  labs(title = "Multivariate Cox regression", x = "Hazard ratio (95% CI)", y = NULL) +
  annotate("text", x = 3.75, y = 4.55, label = "HR (95% CI)", hjust = 0, fontface = "bold", size = 3.7) +
  annotate("text", x = 5.4, y = 4.55, label = "P value", hjust = 0, fontface = "bold", size = 3.7) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.text.y = element_text(size = 11), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 12), plot.margin = margin(10, 25, 10, 10))
ggsave("Fig5D_ENO1_Cox_forest_vector_clean.pdf", p, device = cairo_pdf, width = 8, height = 4.6, units = "in")
cat("Saved figure and check file.\n")
