# ============================================================
# Revision Analyses — JTM Major Revision
# ============================================================

library(Seurat)
library(AUCell)
library(spacexr)
library(ggplot2)
library(dplyr)
library(tidyr)

load("D:/scRNA_project/workspace_backup.RData")

glycolysis_genes <- c("HK1","HK2","GPI","PFKL","PFKP","PFKM",
                      "ALDOA","ALDOB","ALDOC","TPI1","GAPDH",
                      "PGK1","PGAM1","ENO1","ENO2","PKM","LDHA",
                      "LDHB","SLC2A1","SLC2A3","PFKFB3","GCK")

# ============================================================
# Issue 2: Leave-One-Out (LOO) Circularity Analysis
# ============================================================

expr_matrix <- GetAssayData(tumor_hep, layer = "counts")
rankings    <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)

loo_cor <- function(exclude_gene, rankings, expr_matrix) {
  gene_set   <- setdiff(glycolysis_genes, exclude_gene)
  auc        <- AUCell_calcAUC(list(loo = gene_set), rankings,
                               aucMaxRank = ceiling(0.05 * nrow(rankings)))
  auc_scores <- as.numeric(getAUC(auc)["loo", ])
  gene_expr  <- as.numeric(expr_matrix[exclude_gene, ])
  cor(gene_expr, auc_scores, method = "spearman")
}

loo_ENO1  <- loo_cor("ENO1",  rankings, expr_matrix)
loo_LDHA  <- loo_cor("LDHA",  rankings, expr_matrix)
loo_GAPDH <- loo_cor("GAPDH", rankings, expr_matrix)

cat("ENO1  LOO rho:", round(loo_ENO1,  3), "\n")  # 0.254
cat("LDHA  LOO rho:", round(loo_LDHA,  3), "\n")  # 0.346
cat("GAPDH LOO rho:", round(loo_GAPDH, 3), "\n")  # 0.374

# Supplementary Figure S20
df_loo <- data.frame(
  Gene  = rep(c("ENO1","LDHA","GAPDH"), each = 2),
  Type  = rep(c("Standard rho","LOO rho"), 3),
  Value = c(0.554, 0.254, 0.537, 0.346, 0.522, 0.374)
)

p_s20 <- ggplot(df_loo, aes(x = Gene, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Standard rho" = "#E64B35", "LOO rho" = "#4DBBD5")) +
  geom_text(aes(label = round(Value, 3)), position = position_dodge(0.6),
            vjust = -0.5, size = 3.5) +
  ylim(0, 0.65) +
  theme_bw(base_size = 12) +
  labs(title = "LOO Spearman Correlation Analysis",
       y = "Spearman rho", x = "", fill = "") +
  theme(legend.position = "bottom")

ggsave("output/FigS20_LOO_comparison.pdf", p_s20, width = 6, height = 5)
ggsave("output/FigS20_LOO_comparison.png", p_s20, width = 6, height = 5, dpi = 300)

# ============================================================
# Issue 5: Spatial Transcriptomics — RCTD Deconvolution
# ============================================================

samples <- c("HCC1R","HCC2R","HCC3R","HCC4R")

# 加载Visium数据
vis_list2 <- list()
for(s in samples) {
  path <- paste0("D:/scRNA_project/GSE238264/", s, "/", s)
  vis_list2[[s]] <- Load10X_Spatial(path, slice = s)
  vis_list2[[s]] <- NormalizeData(vis_list2[[s]])
  vis_list2[[s]] <- AddModuleScore(vis_list2[[s]],
                                   features = list(glycolysis_genes),
                                   name = "Glycolysis", nbin = 10)
  cat(s, "loaded\n")
}

# 构建scRNA参考（注意替换celltype中的"/"）
counts_ref <- GetAssayData(seu, assay = "RNA", layer = "counts")
cell_types  <- as.factor(gsub("/", "_", seu$celltype))
names(cell_types) <- colnames(seu)
reference   <- Reference(counts_ref, cell_types)
# 注：spacexr会自动下采样每类细胞到10000，属正常行为

# 对四个样本跑RCTD
rctd_results <- list()
for(s in samples) {
  cat("Running RCTD for", s, "...\n")
  obj       <- vis_list2[[s]]
  counts_sp <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  coords    <- GetTissueCoordinates(obj)[, c("x","y")]
  sp_obj    <- SpatialRNA(coords, counts_sp)
  rctd      <- create.RCTD(sp_obj, reference, max_cores = 4)
  rctd      <- run.RCTD(rctd, doublet_mode = "full")
  rctd_results[[s]] <- rctd
  cat(s, "done\n")
}

# 提取结果并比较细胞类型比例
results_list <- list()
for(s in samples) {
  rctd         <- rctd_results[[s]]
  weights      <- as.data.frame(rctd@results$weights)
  weights_norm <- weights / rowSums(weights)
  
  obj          <- vis_list2[[s]]
  glyco_score  <- obj$Glycolysis1
  glyco_group  <- ifelse(glyco_score > median(glyco_score), "GlycoHigh","GlycoLow")
  
  weights_norm$group  <- glyco_group[rownames(weights_norm)]
  weights_norm$sample <- s
  results_list[[s]]   <- weights_norm
}

all_results <- bind_rows(results_list)

# Wilcoxon检验：Hepatocyte比例
hep_high <- all_results %>% filter(group == "GlycoHigh") %>% pull(Hepatocyte)
hep_low  <- all_results %>% filter(group == "GlycoLow")  %>% pull(Hepatocyte)
wt       <- wilcox.test(hep_high, hep_low)
cat("Hepatocyte GlycoHigh mean:", round(mean(hep_high), 3),  # 0.629
    "GlycoLow mean:", round(mean(hep_low), 3),               # 0.670
    "p =", format(wt$p.value, scientific = TRUE), "\n")      # 2.75e-18

# Supplementary Figure S21
plot_df <- all_results %>%
  select(sample, group, Hepatocyte, Myeloid, T_NK, Fibroblast, Endothelial, B) %>%
  pivot_longer(cols = -c(sample, group), names_to = "CellType", values_to = "Proportion")

p_s21 <- ggplot(plot_df, aes(x = group, y = Proportion, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~CellType, scales = "free_y") +
  scale_fill_manual(values = c("GlycoHigh" = "#E64B35","GlycoLow" = "#4DBBD5")) +
  theme_bw(base_size = 12) +
  labs(title = "RCTD Cell Type Composition: GlycoHigh vs GlycoLow Spots",
       x = "", y = "Estimated Proportion") +
  theme(legend.position = "bottom")

ggsave("output/FigS21_RCTD_composition.pdf", p_s21, width = 10, height = 8)
ggsave("output/FigS21_RCTD_composition.png", p_s21, width = 10, height = 8, dpi = 300)
