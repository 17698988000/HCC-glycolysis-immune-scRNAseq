# ============================================================
# 09_TF_activity.R
# Purpose: Quantify transcription factor regulon activity
#          using AUCell in tumor hepatocytes
# Corresponds to: Methods Section 2.12, Results Section 3.11
# ============================================================

library(Seurat)
library(AUCell)
library(ggplot2)

# ── Load data ──────────────────────────────────────────────
seu <- readRDS("D:/scRNA_project/seurat_final.rds")

tumor_samples <- grep("T$", unique(seu$orig.ident), value = TRUE)
seu_hep <- subset(seu, cell_type == "Hepatocyte" & orig.ident %in% tumor_samples)
cat("Tumor hepatocytes:", ncol(seu_hep), "\n")

counts_hep <- GetAssayData(seu_hep, assay = "RNA", layer = "counts")

# ── Define TF target gene sets (from ENCODE ChIP-seq) ─────
tf_targets <- list(
  HIF1A = c("ENO1","LDHA","SLC2A1","PFKFB3","ALDOA","PKM",
             "PGK1","VEGFA","HK2","PFKL","BNIP3","PDK1"),
  MYC   = c("ENO1","LDHA","GPI","TPI1","PFKM","PKM",
             "GAPDH","SLC2A3","PFKL","ALDOA","CDK4","MCM7"),
  SP1   = c("ENO1","HK1","PFKL","GPI","PGAM1","LDHA",
             "SLC2A1","ENO2","PKM","TPI1","VEGFA","CCND1"),
  YAP1  = c("CYR61","CTGF","ANKRD1","CCN2","AMOTL2",
             "BIRC5","AREG","FGF1","MYC","TEAD1","LATS2","NF2"),
  STAT3 = c("MYC","BCL2","VEGFA","HIF1A","SPP1","MMP9",
             "TWIST1","CCND1","IL6R","SOCS3","MCL1","SURVIVIN")
)

# ── AUCell scoring ─────────────────────────────────────────
cat("Building AUCell rankings...\n")
rankings <- AUCell_buildRankings(counts_hep, nCores = 1, plotStats = FALSE)

cat("Calculating TF regulon activity scores...\n")
auc_tf <- AUCell_calcAUC(
  tf_targets, rankings,
  aucMaxRank = ceiling(0.05 * nrow(counts_hep))
)

auc_df <- as.data.frame(t(getAUC(auc_tf)))

# ── Glycolysis AUCell scores ───────────────────────────────
glyco_genes <- list(Glycolysis = c(
  "HK1","HK2","GPI","PFKL","PFKP","PFKM","ALDOA","TPI1","GAPDH",
  "PGK1","PGAM1","ENO1","ENO2","PKM","LDHA","SLC2A1","SLC2A3","PFKFB3"
))
auc_glyco <- AUCell_calcAUC(
  glyco_genes, rankings,
  aucMaxRank = ceiling(0.05 * nrow(counts_hep))
)
auc_df$glycolysis <- as.numeric(getAUC(auc_glyco)["Glycolysis", ])
auc_df$ENO1       <- as.numeric(GetAssayData(seu_hep, layer = "data")["ENO1", ])

# ── Spearman correlations ──────────────────────────────────
tfs <- c("HIF1A","MYC","SP1","YAP1","STAT3")

cor_glyco <- sapply(tfs, function(tf)
  cor(auc_df[[tf]], auc_df$glycolysis, method = "spearman"))
cor_eno1  <- sapply(tfs, function(tf)
  cor(auc_df[[tf]], auc_df$ENO1, method = "spearman"))

cat("\nSpearman rho with glycolysis score:\n"); print(round(cor_glyco, 3))
cat("Spearman rho with ENO1 expression:\n");  print(round(cor_eno1, 3))

# ── Plot FigS16 ────────────────────────────────────────────
plot_df <- data.frame(
  TF          = rep(tfs, 2),
  Correlation = c(cor_glyco, cor_eno1),
  Target      = rep(c("Glycolysis score","ENO1 expression"), each = 5)
)
plot_df$TF <- factor(plot_df$TF, levels = c("HIF1A","MYC","SP1","STAT3","YAP1"))

p <- ggplot(plot_df, aes(x = TF, y = Correlation, fill = Target)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("Glycolysis score" = "#E64B35",
                                "ENO1 expression" = "#4DBBD5")) +
  scale_y_continuous(limits = c(-0.1, 1.0), breaks = seq(0, 1, 0.2)) +
  labs(x = NULL, y = "Spearman correlation (rho)",
       title = "TF Regulon Activity vs. Glycolytic State",
       subtitle = paste0("AUCell-based scoring | Tumor hepatocytes (n = ",
                         ncol(seu_hep), ")"),
       fill = NULL) +
  theme_classic(base_size = 13) +
  theme(legend.position = c(0.75, 0.9),
        axis.text.x = element_text(size = 12, face = "bold"))

ggsave("D:/scRNA_project/FigS16_TF_activity_v2.png",
       p, width = 8, height = 6, dpi = 300)
ggsave("D:/scRNA_project/FigS16_TF_activity_v2.pdf",
       p, width = 8, height = 6)

cat("TF activity analysis complete.\n")
