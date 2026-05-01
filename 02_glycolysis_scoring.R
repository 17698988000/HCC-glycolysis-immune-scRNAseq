# ============================================================
# Script 02: Glycolysis Activity Scoring (AUCell)
#            ENO1-Glycolysis Correlation
# ============================================================

library(Seurat)
library(AUCell)
library(tidyverse)
library(ggpubr)

setwd("D:/scRNA_project")

seu <- readRDS("seurat_8patients.rds")

# --- Glycolysis gene set (22 core genes from HALLMARK_GLYCOLYSIS) ---
glycolysis_genes <- c(
  "HK1","HK2","GPI","PFKL","PFKP","PFKM",
  "ALDOA","ALDOB","ALDOC","TPI1","GAPDH",
  "PGK1","PGAM1","ENO1","ENO2","PKM","LDHA",
  "LDHB","SLC2A1","SLC2A3","PFKFB3","GCK"
)

# --- AUCell scoring ---
expr_matrix <- GetAssayData(seu, layer = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix, 
                                       nCores = 1, 
                                       plotStats = FALSE)
geneSets <- list(glycolysis = glycolysis_genes)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

seu$glycolysis_score <- as.numeric(getAUC(cells_AUC)["glycolysis", ])

# --- Subset tumor hepatocytes and stratify ---
tumor_hepa <- subset(seu, 
                     subset = cell_type == "Hepatocytes" & site == "Tumor")

median_score <- median(tumor_hepa$glycolysis_score)
tumor_hepa$glyco_group <- ifelse(
  tumor_hepa$glycolysis_score >= median_score,
  "TumorGlycoHigh", "TumorGlycoLow"
)
cat("Median glycolysis score:", round(median_score, 4), "\n")
cat("TumorGlycoHigh:", sum(tumor_hepa$glyco_group == "TumorGlycoHigh"), "\n")
cat("TumorGlycoLow:",  sum(tumor_hepa$glyco_group == "TumorGlycoLow"),  "\n")

# --- ENO1-glycolysis correlation (all cells) ---
eno1_all <- as.numeric(GetAssayData(seu, layer = "data")["ENO1", ])
cor_all <- cor.test(eno1_all, seu$glycolysis_score, method = "pearson")
cat("\nENO1 ~ Glycolysis (all cells): R =", 
    round(cor_all$estimate, 3), ", p =", cor_all$p.value, "\n")

# --- ENO1-glycolysis correlation (tumor hepatocytes only) ---
eno1_hepa <- as.numeric(
  GetAssayData(tumor_hepa, layer = "data")["ENO1", ]
)
cor_hepa <- cor.test(eno1_hepa, tumor_hepa$glycolysis_score, 
                     method = "pearson")
cat("ENO1 ~ Glycolysis (tumor hepatocytes): R =",
    round(cor_hepa$estimate, 3), ", p =", cor_hepa$p.value, "\n")

# --- Sensitivity: ENO1-excluded AUCell ---
glycolysis_no_eno1 <- setdiff(glycolysis_genes, "ENO1")
geneSets_noENO1 <- list(glycolysis_noENO1 = glycolysis_no_eno1)
cells_AUC_noENO1 <- AUCell_calcAUC(geneSets_noENO1, cells_rankings,
                                   aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))
score_noENO1 <- as.numeric(getAUC(cells_AUC_noENO1)["glycolysis_noENO1",
                                                    colnames(tumor_hepa)])
cor_sensitivity <- cor.test(eno1_hepa, score_noENO1, method = "pearson")
cat("ENO1 ~ ENO1-excluded score: R =",
    round(cor_sensitivity$estimate, 3), "\n")

# --- Save updated object ---
saveRDS(seu, "seurat_final.rds")
saveRDS(tumor_hepa, "tumor_hepatocytes.rds")
cat("Saved: seurat_final.rds, tumor_hepatocytes.rds\n")