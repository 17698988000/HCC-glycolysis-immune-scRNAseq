# ============================================================
# Script 04: TCGA-LIHC Bulk Validation
#            - ENO1 expression (tumor vs normal)
#            - Survival analysis
#            - Immune infiltration correlation
# ============================================================

library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(org.Hs.eg.db)

setwd("D:/scRNA_project")

# --- Load TCGA data ---
fpkm <- read_tsv("TCGA-LIHC.star_fpkm.tsv.gz")
gene_names <- fpkm[[1]]
expr_mat <- as.matrix(fpkm[, -1])
rownames(expr_mat) <- gene_names

# Convert Ensembl to gene symbols
ensembl_ids <- gsub("\\..*", "", rownames(expr_mat))
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = ensembl_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
rownames(expr_mat) <- ifelse(is.na(gene_symbols),
                             ensembl_ids, gene_symbols)

# Separate tumor vs normal
tumor_samples  <- colnames(expr_mat)[substr(colnames(expr_mat),14,15)=="01"]
normal_samples <- colnames(expr_mat)[substr(colnames(expr_mat),14,15)=="11"]
expr_tumor  <- expr_mat[, tumor_samples]
expr_normal <- expr_mat[, normal_samples]

expr_log2 <- log2(expr_tumor + 1)

# --- ENO1 differential expression ---
eno1_tumor  <- log2(as.numeric(expr_mat["ENO1", tumor_samples])  + 1)
eno1_normal <- log2(as.numeric(expr_mat["ENO1", normal_samples]) + 1)
wilcox_result <- wilcox.test(eno1_tumor, eno1_normal)
cat("ENO1 tumor vs normal, p =", wilcox_result$p.value, "\n")

# --- Survival analysis ---
survival_data <- read_tsv("TCGA-LIHC.survival.tsv.gz")
clinical_data <- read_tsv("TCGA-LIHC.clinical.tsv.gz")

surv_df <- data.frame(
  sample    = names(eno1_tumor),
  ENO1_expr = eno1_tumor
) %>%
  left_join(survival_data, by = c("sample" = "sample")) %>%
  drop_na(OS, OS.time)

# Optimal cutpoint
cut_result <- surv_cutpoint(surv_df, time = "OS.time",
                            event = "OS", variables = "ENO1_expr")
optimal_cut <- summary(cut_result)$cutpoint
surv_df$ENO1_group <- ifelse(surv_df$ENO1_expr >= optimal_cut,
                             "High", "Low")

km_fit <- survfit(Surv(OS.time, OS) ~ ENO1_group, data = surv_df)
logrank <- survdiff(Surv(OS.time, OS) ~ ENO1_group, data = surv_df)
cat("KM log-rank p =",
    1 - pchisq(logrank$chisq, df = 1), "\n")

# --- Immune infiltration correlation ---
sample_names <- colnames(expr_log2)
eno1_vec <- as.numeric(expr_log2["ENO1", ])
names(eno1_vec) <- sample_names

immune_markers <- list(
  CD8T      = c("CD8A","CD8B","GZMB","PRF1","IFNG"),
  CD4T      = c("CD4","IL7R","TCF7"),
  NK        = c("KLRD1","GNLY","NKG7","KLRB1"),
  M1_Mac    = c("CD86","IL1B","TNF","NOS2"),
  M2_Mac    = c("CD163","MRC1","CCL18","TGFB1"),
  MDSC      = c("S100A9","S100A8","LYZ","VCAN"),
  Treg      = c("FOXP3","IL2RA","CTLA4"),
  Exhausted = c("PDCD1","HAVCR2","TIGIT","LAG3","TOX")
)

expr_zscore <- t(scale(t(expr_log2)))
colnames(expr_zscore) <- sample_names

immune_genes_all <- unique(unlist(immune_markers))
immune_genes_present <- intersect(immune_genes_all, rownames(expr_log2))

cor_results <- lapply(immune_genes_present, function(gene) {
  test <- cor.test(eno1_vec,
                   as.numeric(expr_log2[gene, ]),
                   method = "spearman")
  data.frame(gene = gene,
             rho  = test$estimate,
             pval = test$p.value)
}) %>% bind_rows()

cor_results$fdr <- p.adjust(cor_results$pval, method = "BH")

cat("\nImmune correlation summary:\n")
cell_label <- stack(immune_markers)
colnames(cell_label) <- c("gene", "cell_type")
cor_results <- left_join(cor_results, cell_label, by = "gene")
cor_results %>%
  group_by(cell_type) %>%
  summarise(mean_rho = round(mean(rho), 3),
            sig_genes = sum(fdr < 0.05)) %>%
  arrange(mean_rho) %>%
  print()

cat("TCGA validation complete.\n")