# ============================================================
# 08_glycolysis_gradient.R
# Purpose: Analyze immunosuppressive ligand expression along
#          the continuous glycolysis activity gradient
# Corresponds to: Methods Section 2.11, Results Section 3.10
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)

# ── Load data ──────────────────────────────────────────────
seu <- readRDS("D:/scRNA_project/seurat_final.rds")

# Extract tumor hepatocytes
tumor_samples <- grep("T$", unique(seu$orig.ident), value = TRUE)
seu_hep <- subset(seu, cell_type == "Hepatocyte" & orig.ident %in% tumor_samples)
cat("Tumor hepatocytes:", ncol(seu_hep), "\n")

# ── Get glycolysis AUCell scores ───────────────────────────
# These should already be in metadata from 02_glycolysis_scoring.R
# If not, re-run AUCell scoring here
if (!"glycolysis_AUC" %in% colnames(seu_hep@meta.data)) {
  library(AUCell)
  glyco_genes <- list(Glycolysis = c(
    "HK1","HK2","GPI","PFKL","PFKP","PFKM","ALDOA","TPI1","GAPDH",
    "PGK1","PGAM1","ENO1","ENO2","PKM","LDHA","SLC2A1","SLC2A3","PFKFB3"
  ))
  counts_hep <- GetAssayData(seu_hep, assay = "RNA", layer = "counts")
  rankings   <- AUCell_buildRankings(counts_hep, nCores = 1, plotStats = FALSE)
  auc_scores <- AUCell_calcAUC(glyco_genes, rankings,
                                aucMaxRank = ceiling(0.05 * nrow(counts_hep)))
  seu_hep$glycolysis_AUC <- as.numeric(getAUC(auc_scores)["Glycolysis", ])
}

# ── 修复：用 ntile() 强制等人数分箱，避免重复分位数导致空 bin ──────
library(dplyr)
library(tidyr)

seu_hep$glyco_bin <- ntile(seu_hep$glycolysis_AUC, 10)
cat("Bin 分布（每个 bin 细胞数）:\n")
print(table(seu_hep$glyco_bin))   # 确认 10 个 bin 都有细胞

# ── 计算每个 bin 的均值 ──────────────────────────────────────────
target_genes <- c("ENO1","LDHA","SPP1","MIF","PTGES")
expr_data    <- GetAssayData(seu_hep, assay = "RNA", layer = "data")

bin_means <- sapply(1:10, function(b) {
  cells <- colnames(seu_hep)[seu_hep$glyco_bin == b]
  rowMeans(expr_data[target_genes, cells, drop = FALSE])
})

gradient_df        <- as.data.frame(t(bin_means))
gradient_df$bin    <- 1:10
rownames(gradient_df) <- NULL
print(gradient_df)   # 确认 10 行都有数据，无 NaN

# ── 出图（Supplementary Figure S15）─────────────────────────────
plot_df <- pivot_longer(gradient_df, cols = all_of(target_genes),
                        names_to = "gene", values_to = "expression")
plot_df$gene_type <- ifelse(plot_df$gene %in% c("ENO1","LDHA"),
                             "Glycolytic enzyme", "Immunosuppressive ligand")

p <- ggplot(plot_df, aes(x = bin, y = expression,
                          color = gene, linetype = gene_type)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:10, labels = paste0("Bin ", 1:10)) +
  scale_linetype_manual(values = c("Glycolytic enzyme" = "solid",
                                   "Immunosuppressive ligand" = "dashed")) +
  labs(x = "Glycolysis activity bin (low → high)",
       y = "Mean log-normalized expression",
       title = "Immunosuppressive Ligand Expression Along Glycolysis Gradient",
       color = "Gene", linetype = "Category") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("D:/scRNA_project/FigS15_glycolysis_gradient.png",
       p, width = 9, height = 5, dpi = 300)
ggsave("D:/scRNA_project/FigS15_glycolysis_gradient.pdf",
       p, width = 9, height = 5)
cat("完成。x 轴应连续显示 Bin 1 到 Bin 10。\n")
