# Figure 5C: ENO1 expression across AJCC pathologic stage
# Output:
#   Fig5C_ENO1_AJCC_stage_vector_clean.pdf
#   Fig5C_ENO1_AJCC_stage_check.txt

source("00_config_paths.R")
ensure_packages(c("readr", "ggplot2", "AnnotationDbi", "org.Hs.eg.db"))
suppressPackageStartupMessages({library(readr); library(ggplot2); library(AnnotationDbi); library(org.Hs.eg.db)})
expr_file <- "TCGA-LIHC.star_fpkm.tsv.gz"
clin_file <- "TCGA-LIHC.clinical.tsv.gz"

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
expr_df <- data.frame(sample16 = substr(tumor_samples, 1, 16),
                      ENO1_expr = as.numeric(log2(expr_mat["ENO1", tumor_samples] + 1)), stringsAsFactors = FALSE)
expr_df <- expr_df[!duplicated(expr_df$sample16), ]

clin_raw <- as.data.frame(readr::read_tsv(clin_file, show_col_types = FALSE))
stage_candidates <- colnames(clin_raw)[grepl("ajcc.*pathologic.*stage|pathologic.*stage|tumor.*stage|stage", colnames(clin_raw), ignore.case = TRUE)]
stage_col <- stage_candidates[grepl("ajcc.*pathologic.*stage", stage_candidates, ignore.case = TRUE)][1]
if (is.na(stage_col)) stage_col <- stage_candidates[grepl("pathologic.*stage", stage_candidates, ignore.case = TRUE)][1]
if (is.na(stage_col)) stage_col <- stage_candidates[1]
if (is.na(stage_col)) stop("No stage column found.")
clin_df <- data.frame(sample16 = substr(as.character(clin_raw$sample), 1, 16), stage_raw = as.character(clin_raw[[stage_col]]), stringsAsFactors = FALSE)
clin_df <- clin_df[!duplicated(clin_df$sample16), ]
plot_df <- merge(expr_df, clin_df, by = "sample16")
stage_upper <- toupper(plot_df$stage_raw)
plot_df$AJCC_stage <- ifelse(grepl("IV", stage_upper), "Stage IV",
                             ifelse(grepl("III", stage_upper), "Stage III",
                                    ifelse(grepl("II", stage_upper), "Stage II",
                                           ifelse(grepl("I", stage_upper), "Stage I", NA))))
plot_df <- plot_df[!is.na(plot_df$AJCC_stage), ]
plot_df$AJCC_stage <- factor(plot_df$AJCC_stage, levels = c("Stage I", "Stage II", "Stage III", "Stage IV"))
kw <- kruskal.test(ENO1_expr ~ AJCC_stage, data = plot_df)
kw_label <- ifelse(kw$p.value < 0.001, "Kruskal-Wallis p < 0.001", paste0("Kruskal-Wallis p = ", formatC(kw$p.value, digits = 3, format = "e")))

check_file <- "Fig5C_ENO1_AJCC_stage_check.txt"
sink(check_file)
cat("===== Fig5C ENO1 AJCC stage check =====\n\n")
cat("Using stage column:", stage_col, "\n")
cat("Matched samples with stage:", nrow(plot_df), "\n")
cat("Collapsed AJCC stage counts:\n"); print(table(plot_df$AJCC_stage))
cat("Kruskal-Wallis p:", kw$p.value, "\n")
sink()

p <- ggplot(plot_df, aes(x = AJCC_stage, y = ENO1_expr)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, color = "grey20", fill = "white") +
  geom_jitter(width = 0.18, size = 1.4, alpha = 0.65) +
  annotate("text", x = 3.0, y = max(plot_df$ENO1_expr, na.rm = TRUE) * 0.98, label = kw_label, size = 4.2) +
  labs(x = "AJCC pathologic stage", y = "ENO1 expression (log2 FPKM + 1)") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(size = 11), axis.title = element_text(size = 12), plot.margin = margin(10, 12, 10, 10))
ggsave("Fig5C_ENO1_AJCC_stage_vector_clean.pdf", p, device = cairo_pdf, width = 5.6, height = 4.8, units = "in")
cat("Saved figure and check file.\n")
