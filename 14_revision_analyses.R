# ============================================================
# OXPHOS代谢特异性偏相关分析（Revision新增）
# ============================================================

OXPHOS_genes <- c(
  "NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
  "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFS1","NDUFS2","NDUFS3",
  "SDHA","SDHB","SDHC","SDHD",
  "UQCRB","UQCRC1","UQCRC2","UQCRQ","CYC1",
  "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX7A2","COX8A",
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
  "ATP5MC1","ATP5MC2","ATP5MC3","ATP5PB"
)

tumor_hep <- subset(seu, site == "Tumor" & celltype == "Hepatocyte")
expr_matrix <- GetAssayData(tumor_hep, assay = "RNA", layer = "counts")
rankings_oxphos <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
auc_oxphos <- AUCell_calcAUC(list(OXPHOS = OXPHOS_genes), rankings_oxphos,
                              aucMaxRank = ceiling(0.05 * nrow(rankings_oxphos)))
tumor_hep$OXPHOS_AUC <- as.numeric(getAUC(auc_oxphos)["OXPHOS", ])

spp1 <- GetAssayData(tumor_hep, layer = "data")["SPP1", ]
mif  <- GetAssayData(tumor_hep, layer = "data")["MIF", ]

library(ppcor)
df_partial <- data.frame(
  SPP1 = as.numeric(spp1), MIF = as.numeric(mif),
  Glycolysis = tumor_hep$Glycolysis_AUC,
  OXPHOS = tumor_hep$OXPHOS_AUC
)

pc_spp1_glyco  <- pcor.test(df_partial$SPP1, df_partial$Glycolysis, df_partial$OXPHOS, method = "spearman")
pc_spp1_oxphos <- pcor.test(df_partial$SPP1, df_partial$OXPHOS, df_partial$Glycolysis, method = "spearman")
pc_mif_glyco   <- pcor.test(df_partial$MIF,  df_partial$Glycolysis, df_partial$OXPHOS, method = "spearman")
pc_mif_oxphos  <- pcor.test(df_partial$MIF,  df_partial$OXPHOS, df_partial$Glycolysis, method = "spearman")

# ============================================================
# DoRothEA全基因组TF活性分析（Revision新增，替换5-TF候选分析）
# ============================================================

library(dorothea)
library(viper)

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B"))

expr_mat <- as.matrix(GetAssayData(tumor_hep, assay = "RNA", layer = "data"))
reg_list <- split(regulons, regulons$tf)
viper_regulons <- lapply(reg_list, function(x) {
  list(tfmode = setNames(x$mor, x$target),
       likelihood = rep(1, nrow(x)))
})

tf_act <- viper(expr_mat, regulon = viper_regulons,
                eset.filter = FALSE, minsize = 4, verbose = FALSE)

glyco_scores <- tumor_hep$Glycolysis_AUC
cor_results <- data.frame(TF = rownames(tf_act), rho = NA, p = NA)
for(i in 1:nrow(tf_act)) {
  ct <- cor.test(tf_act[i, ], glyco_scores, method = "spearman")
  cor_results$rho[i] <- ct$estimate
  cor_results$p[i]   <- ct$p.value
}
cor_results$padj <- p.adjust(cor_results$p, method = "BH")
cor_results <- cor_results[order(-abs(cor_results$rho)), ]
