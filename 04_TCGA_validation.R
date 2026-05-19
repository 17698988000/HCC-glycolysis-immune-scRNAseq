# ============================================================
# Script 04: TCGA-LIHC Bulk Validation
# Sections 2.7 & 2.8
# - ENO1 tumor vs normal expression
# - Survival analysis (primary: median cutoff; sensitivity: optimal cutoff)
# - Multivariate Cox regression
# - Immune infiltration correlation (panel mean method)
# - ssGSEA (22 immune cell types, Supplementary Figure S11)
# - TIDE T cell dysfunction / exclusion scores (Supplementary Figure S13)
# ============================================================

library(tidyverse)
library(survival)
library(survminer)
library(ggpubr)
library(org.Hs.eg.db)
library(GSVA)          # for ssGSEA

setwd("D:/scRNA_project")

# ============================================================
# 1. Load and preprocess TCGA-LIHC data
# ============================================================
cat("=== Loading TCGA-LIHC data ===\n")

fpkm      <- read_tsv("TCGA-LIHC.star_fpkm.tsv.gz")
gene_names <- fpkm[[1]]
expr_mat   <- as.matrix(fpkm[, -1])
rownames(expr_mat) <- gene_names

# Convert Ensembl IDs to gene symbols
ensembl_ids  <- gsub("\\..*", "", rownames(expr_mat))
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys     = ensembl_ids,
                       column   = "SYMBOL",
                       keytype  = "ENSEMBL",
                       multiVals = "first")
rownames(expr_mat) <- ifelse(is.na(gene_symbols), ensembl_ids, gene_symbols)

# Separate tumor (01) vs normal (11) samples
tumor_samples  <- colnames(expr_mat)[substr(colnames(expr_mat), 14, 15) == "01"]
normal_samples <- colnames(expr_mat)[substr(colnames(expr_mat), 14, 15) == "11"]
cat("Tumor samples:", length(tumor_samples),
    "| Normal samples:", length(normal_samples), "\n")

expr_tumor  <- expr_mat[, tumor_samples]
expr_normal <- expr_mat[, normal_samples]

# Log2 transform
expr_log2_tumor  <- log2(expr_tumor  + 1)
expr_log2_normal <- log2(expr_normal + 1)

# ============================================================
# 2. ENO1 differential expression: tumor vs normal
# ============================================================
cat("\n=== ENO1: tumor vs normal ===\n")

eno1_tumor  <- as.numeric(expr_log2_tumor ["ENO1", ])
eno1_normal <- as.numeric(expr_log2_normal["ENO1", ])

wt <- wilcox.test(eno1_tumor, eno1_normal)
cat("ENO1 tumor mean:", round(mean(eno1_tumor), 3),
    "| normal mean:", round(mean(eno1_normal), 3),
    "| Wilcoxon p =", formatC(wt$p.value, format = "e", digits = 2), "\n")

# ============================================================
# 3. Survival analysis
# ============================================================
cat("\n=== Survival analysis ===\n")

survival_data <- read_tsv("TCGA-LIHC.survival.tsv.gz")
clinical_data <- read_tsv("TCGA-LIHC.clinical.tsv.gz")

surv_df <- data.frame(
  sample    = tumor_samples,
  ENO1_expr = eno1_tumor
) %>%
  left_join(survival_data, by = c("sample" = "sample")) %>%
  left_join(clinical_data %>%
              select(case_submitter_id, age_at_index, gender,
                     ajcc_pathologic_stage) %>%
              mutate(sample = paste0(case_submitter_id, "-01A")),
            by = "sample") %>%
  drop_na(OS, OS.time)

cat("Samples with survival data:", nrow(surv_df), "\n")

# ── FIX: Primary analysis uses MEDIAN cutoff ──
median_cut          <- median(surv_df$ENO1_expr)
surv_df$ENO1_group  <- ifelse(surv_df$ENO1_expr >= median_cut, "High", "Low")
cat("Median cutoff:", round(median_cut, 3),
    "| High:", sum(surv_df$ENO1_group == "High"),
    "| Low:", sum(surv_df$ENO1_group == "Low"), "\n")

# Primary KM (median cutoff) — main result, Supplementary Figure S4
km_fit_median  <- survfit(Surv(OS.time, OS) ~ ENO1_group, data = surv_df)
lr_median      <- survdiff(Surv(OS.time, OS) ~ ENO1_group, data = surv_df)
pval_median    <- 1 - pchisq(lr_median$chisq, df = 1)
cat("Primary KM (median cutoff) log-rank p =",
    formatC(pval_median, format = "e", digits = 2), "\n")   # expect 3.31e-5

p_km_primary <- ggsurvplot(
  km_fit_median, data = surv_df,
  pval = TRUE, risk.table = TRUE,
  palette = c("#E64B35", "#4DBBD5"),
  legend.labs = c("ENO1-High", "ENO1-Low"),
  xlab = "Time (days)", ylab = "Overall survival probability",
  title = "ENO1 Overall Survival (Primary: median cutoff)",
  ggtheme = theme_classic(base_size = 13)
)
ggsave("FigS4_ENO1_KM_median_cutoff.pdf",
       print(p_km_primary), width = 8, height = 7)

# Sensitivity KM (optimal cutoff) — Figure 5B inset
cut_result   <- surv_cutpoint(surv_df, time = "OS.time",
                               event = "OS", variables = "ENO1_expr",
                               minprop = 0.1)
optimal_cut  <- summary(cut_result)$cutpoint
surv_df$ENO1_group_opt <- ifelse(surv_df$ENO1_expr >= optimal_cut, "High", "Low")
cat("Optimal cutoff:", round(optimal_cut, 3),
    "| High:", sum(surv_df$ENO1_group_opt == "High"),
    "| Low:", sum(surv_df$ENO1_group_opt == "Low"), "\n")

km_fit_opt <- survfit(Surv(OS.time, OS) ~ ENO1_group_opt, data = surv_df)
lr_opt     <- survdiff(Surv(OS.time, OS) ~ ENO1_group_opt, data = surv_df)
pval_opt   <- 1 - pchisq(lr_opt$chisq, df = 1)
cat("Sensitivity KM (optimal cutoff) log-rank p =",
    formatC(pval_opt, format = "e", digits = 2), "\n")     # expect 0.019

p_km_sens <- ggsurvplot(
  km_fit_opt, data = surv_df,
  pval = TRUE, risk.table = TRUE,
  palette = c("#E64B35", "#4DBBD5"),
  legend.labs = c("ENO1-High", "ENO1-Low"),
  xlab = "Time (days)", ylab = "Overall survival probability",
  title = paste0("ENO1 Survival (Sensitivity: optimal cutoff = ",
                 round(optimal_cut, 3), ")"),
  ggtheme = theme_classic(base_size = 13)
)
ggsave("Fig5B_ENO1_KM_optimal_cutoff.pdf",
       print(p_km_sens), width = 8, height = 7)

# ============================================================
# 4. Multivariate Cox regression for ENO1 (Section 2.7)
# ============================================================
cat("\n=== Multivariate Cox regression (ENO1) ===\n")

# Prepare clinical covariates — filter to samples with complete staging data
surv_cox <- surv_df %>%
  mutate(
    stage_binary = ifelse(
      grepl("III|IV", ajcc_pathologic_stage), "III-IV", "I-II"
    ),
    age_group = ifelse(age_at_index >= 60, ">=60", "<60")
  ) %>%
  drop_na(stage_binary, age_at_index, gender)

cat("Samples for multivariate Cox:", nrow(surv_cox), "\n")  # expect ~374

# Use optimal cutoff for Cox (consistent with Figure 5D in manuscript)
surv_cox$ENO1_group_opt <- ifelse(
  surv_cox$ENO1_expr >= optimal_cut, "High", "Low"
)
surv_cox$ENO1_group_opt <- relevel(factor(surv_cox$ENO1_group_opt), ref = "Low")
surv_cox$stage_binary   <- relevel(factor(surv_cox$stage_binary),   ref = "I-II")

# Univariate Cox
cox_uni <- coxph(Surv(OS.time, OS) ~ ENO1_group_opt, data = surv_cox)
cat("Univariate Cox:\n"); print(summary(cox_uni)$conf.int)

# Multivariate Cox (adjusted for age, gender, AJCC stage)
cox_multi <- coxph(
  Surv(OS.time, OS) ~ ENO1_group_opt + age_at_index + gender + stage_binary,
  data = surv_cox
)
cat("Multivariate Cox:\n"); print(summary(cox_multi)$conf.int)
cat("C-index:", summary(cox_multi)$concordance[1], "\n")   # expect ~0.623

# Forest plot — Figure 5D
png("Fig5D_ENO1_Cox_forest.png", width = 8, height = 5,
    units = "in", res = 300, bg = "white")
ggforest(cox_multi, data = surv_cox,
         main = "Multivariate Cox: ENO1 and Clinical Covariates")
dev.off()

# ============================================================
# 5. Immune infiltration correlation (panel mean method, Section 2.8)
# ============================================================
cat("\n=== Immune infiltration correlation ===\n")

# Filter to 371 tumor samples with complete expression data
eno1_vec    <- as.numeric(expr_log2_tumor["ENO1", ])
names(eno1_vec) <- tumor_samples

immune_panels <- list(
  CD8T      = c("CD8A", "CD8B", "GZMB", "PRF1", "IFNG"),
  CD4T      = c("CD4",  "IL7R", "TCF7"),
  NK        = c("KLRD1","GNLY", "NKG7", "KLRB1"),
  M1_Mac    = c("CD86", "IL1B", "TNF",  "NOS2"),
  M2_Mac    = c("CD163","MRC1", "CCL18","TGFB1"),
  MDSC      = c("S100A9","S100A8","LYZ", "VCAN"),
  Treg      = c("FOXP3","IL2RA","CTLA4"),
  Exhausted = c("PDCD1","HAVCR2","TIGIT","LAG3","TOX")
)

# ── FIX: Compute PANEL MEAN per sample, then correlate with ENO1 ──
cor_results <- lapply(names(immune_panels), function(ct) {
  genes        <- intersect(immune_panels[[ct]], rownames(expr_log2_tumor))
  panel_means  <- colMeans(expr_log2_tumor[genes, , drop = FALSE])
  # Keep only samples present in both vectors
  common       <- intersect(names(eno1_vec), names(panel_means))
  test         <- cor.test(eno1_vec[common], panel_means[common],
                           method = "spearman")
  data.frame(
    cell_type = ct,
    rho       = round(test$estimate, 3),
    pval      = test$p.value,
    n_genes   = length(genes)
  )
}) %>% bind_rows()

cor_results$fdr <- p.adjust(cor_results$pval, method = "BH")
cor_results$sig <- ifelse(cor_results$fdr < 0.05, "*", "ns")
cat("\nImmune panel correlations with ENO1:\n")
print(cor_results)

# Figure 5F
p_immune <- ggplot(cor_results,
                   aes(x = reorder(cell_type, rho), y = rho,
                       fill = rho > 0)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = sig,
                y = rho + sign(rho) * 0.005),
            size = 5, vjust = 0) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                    guide = "none") +
  coord_flip() +
  labs(x = NULL, y = "Spearman rho",
       title = "ENO1 vs. Immune Cell Signatures (TCGA-LIHC)",
       subtitle = "Panel mean expression per sample; * FDR < 0.05") +
  theme_classic(base_size = 12)

ggsave("Fig5F_ENO1_immune_correlation.pdf", p_immune, width = 7, height = 5)

# ============================================================
# 6. ssGSEA — 22 immune cell types (Supplementary Figure S11)
# ============================================================
cat("\n=== ssGSEA immune infiltration ===\n")

# LM22-derived marker gene sets (22 immune cell types)
lm22_sets <- list(
  B_cells_naive         = c("CD19","MS4A1","CD79A","IGHM","IGHD"),
  B_cells_memory        = c("CD19","MS4A1","CD27","IGHG1"),
  Plasma_cells          = c("MZB1","SDC1","IGHG1","IGKC","IGLC2"),
  T_cells_CD8           = c("CD8A","CD8B","GZMB","PRF1"),
  T_cells_CD4_naive     = c("CD4","CCR7","SELL","TCF7"),
  T_cells_CD4_memory_resting = c("CD4","IL7R","S100B"),
  T_cells_CD4_memory_activated = c("CD4","TNFRSF4","ICOS","CD69"),
  T_cells_follicular    = c("CXCR5","BCL6","PDCD1","ICOS"),
  T_cells_regulatory    = c("FOXP3","IL2RA","CTLA4","IKZF2"),
  T_cells_gamma_delta   = c("TRDC","TRGC1","TRGC2","KLRD1"),
  NK_cells_resting      = c("KLRD1","KLRF1","NCAM1","NKG7"),
  NK_cells_activated    = c("KLRD1","GZMB","PRF1","IFNG"),
  Monocytes             = c("CD14","LYZ","S100A9","FCN1"),
  Macrophages_M0        = c("CD68","ADCP","MRC1","CD163"),
  Macrophages_M1        = c("CD86","IL1B","TNF","NOS2","CXCL10"),
  Macrophages_M2        = c("CD163","MRC1","CCL18","TGFB1","IL10"),
  Dendritic_cells_resting   = c("ITGAX","CD1C","CLEC9A","FCER1A"),
  Dendritic_cells_activated = c("ITGAX","CD80","CD86","CCR7"),
  Mast_cells_resting    = c("KIT","MS4A2","TPSAB1","CPA3"),
  Mast_cells_activated  = c("KIT","MS4A2","IL13","TPSAB1"),
  Eosinophils           = c("CCR3","IL5RA","SIGLEC8","EPX"),
  Neutrophils           = c("FCGR3B","CSF3R","CXCR2","FFAR2")
)

# Filter gene sets to those present in expression matrix
lm22_filtered <- lapply(lm22_sets, function(g)
  intersect(g, rownames(expr_mat))
)
lm22_filtered <- lm22_filtered[sapply(lm22_filtered, length) >= 3]
cat("Gene sets with >= 3 genes:", length(lm22_filtered), "\n")

# Run ssGSEA on all 424 samples (374 tumor + 50 normal)
expr_all_log2 <- log2(expr_mat + 1)
ssgsea_scores <- gsva(as.matrix(expr_all_log2),
                      lm22_filtered,
                      method     = "ssgsea",
                      kcdf       = "Gaussian",
                      min.sz     = 3,
                      verbose    = FALSE)

# Correlate with ENO1 across all 424 samples
eno1_all   <- as.numeric(expr_all_log2["ENO1", ])
ssgsea_cor <- apply(ssgsea_scores, 1, function(score) {
  ct <- cor.test(eno1_all, score, method = "spearman")
  c(rho = ct$estimate, pval = ct$p.value)
})
ssgsea_cor_df           <- as.data.frame(t(ssgsea_cor))
ssgsea_cor_df$cell_type <- rownames(ssgsea_cor_df)
ssgsea_cor_df$padj      <- p.adjust(ssgsea_cor_df$pval, method = "BH")

cat("\nssGSEA correlations (significant, FDR<0.05):\n")
print(ssgsea_cor_df[ssgsea_cor_df$padj < 0.05,
                    c("cell_type","rho","padj")])

write.csv(ssgsea_cor_df, "FigS11_ssGSEA_ENO1_correlations.csv",
          row.names = FALSE)

# Supplementary Figure S11
p_s11 <- ggplot(ssgsea_cor_df,
                aes(x = reorder(cell_type, rho), y = rho,
                    fill = padj < 0.05)) +
  geom_bar(stat = "identity", width = 0.75) +
  scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey70"),
                    name = "FDR < 0.05") +
  coord_flip() +
  labs(x = NULL, y = "Spearman rho",
       title = "ENO1 vs. 22 Immune Cell Types (ssGSEA, TCGA-LIHC n=424)",
       subtitle = "Red = FDR < 0.05") +
  theme_classic(base_size = 11)

ggsave("FigS11_ssGSEA_ENO1_immune.pdf",  p_s11, width = 8, height = 7)
ggsave("FigS11_ssGSEA_ENO1_immune.png",  p_s11, width = 8, height = 7, dpi = 300)
cat("Supplementary Figure S11 saved.\n")

# ============================================================
# 7. TIDE T cell dysfunction and exclusion scores
#    (Supplementary Figure S13, Section 2.8)
# ============================================================
cat("\n=== TIDE immune evasion scores ===\n")

# Restrict to 374 tumor samples only (normal tissue excluded)
expr_tumor_log2 <- log2(expr_tumor + 1)

# T cell dysfunction markers (9 genes)
dysfunction_genes <- c("PDCD1","HAVCR2","TIGIT","LAG3","TOX",
                        "CTLA4","CD244","ENTPD1","CD160")

# T cell exclusion markers (13 genes)
exclusion_genes <- c("TGFB1","TGFB2","TGFB3","THBS1","THBS2",
                      "MMP2","MMP9","FAP","FN1","VIM",
                      "CDH2","ZEB1","SNAI1")

# Compute signature means per sample
dysfunc_genes_present <- intersect(dysfunction_genes, rownames(expr_tumor_log2))
exclusion_genes_present <- intersect(exclusion_genes, rownames(expr_tumor_log2))

cat("Dysfunction genes present:", length(dysfunc_genes_present), "/ 9\n")
cat("Exclusion genes present:  ", length(exclusion_genes_present), "/ 13\n")

tide_df <- data.frame(
  sample     = tumor_samples,
  ENO1_expr  = as.numeric(expr_tumor_log2["ENO1", ]),
  Dysfunction = colMeans(expr_tumor_log2[dysfunc_genes_present, ]),
  Exclusion   = colMeans(expr_tumor_log2[exclusion_genes_present, ])
) %>%
  mutate(
    TIDE_score = Exclusion - Dysfunction,
    ENO1_group = ifelse(ENO1_expr >= median(ENO1_expr), "High", "Low")
  )

cat("ENO1 High n:", sum(tide_df$ENO1_group == "High"),
    "| Low n:", sum(tide_df$ENO1_group == "Low"), "\n")

wt_dysfunc  <- wilcox.test(Dysfunction ~ ENO1_group, data = tide_df)
wt_exclusion <- wilcox.test(Exclusion  ~ ENO1_group, data = tide_df)
wt_tide      <- wilcox.test(TIDE_score ~ ENO1_group, data = tide_df)
cat("Dysfunction p =",  formatC(wt_dysfunc$p.value,   format = "e", digits = 2),
    "\n")  # expect 0.0001
cat("Exclusion p =",    formatC(wt_exclusion$p.value, format = "e", digits = 2),
    "\n")  # expect ~0.095 (ns)
cat("TIDE p =",         formatC(wt_tide$p.value,      format = "e", digits = 2),
    "\n")  # expect ~0.017

# Supplementary Figure S13
tide_long <- tide_df %>%
  select(ENO1_group, Dysfunction, Exclusion, TIDE_score) %>%
  pivot_longer(cols = -ENO1_group,
               names_to = "Score", values_to = "Value")

p_s13 <- ggboxplot(
  tide_long,
  x = "ENO1_group", y = "Value", fill = "ENO1_group",
  facet.by = "Score", scales = "free_y",
  palette = c("High" = "#E64B35", "Low" = "#4DBBD5"),
  xlab = "ENO1 Expression Group",
  title = "TIDE Immune Evasion Scores by ENO1 Expression (TCGA-LIHC)"
) +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     comparisons = list(c("High", "Low"))) +
  theme(legend.position = "none")

ggsave("FigS13_TIDE_scores.pdf", p_s13, width = 9, height = 5)
ggsave("FigS13_TIDE_scores.png", p_s13, width = 9, height = 5, dpi = 300)
cat("Supplementary Figure S13 saved.\n")

cat("\n=== Script 04 complete ===\n")

# Final submission clean vector versions of Fig5B-Fig6B are regenerated in:
# revision_figure_restore/plot_Fig5B_ENO1_KM_optimal.R
# revision_figure_restore/plot_Fig5C_ENO1_AJCC_stage.R
# revision_figure_restore/plot_Fig5D_ENO1_Cox_forest.R
# revision_figure_restore/plot_Fig6A_four_gene_score_KM.R
# revision_figure_restore/plot_Fig6B_four_gene_score_Cox_forest.R
