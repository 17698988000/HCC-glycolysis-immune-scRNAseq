# ============================================================
# Script 05: Independent Validation in GSE125449
#            - ENO1-glycolysis correlation
#            - SPP1 and MIF upregulation in GlycoHigh
# ============================================================

library(Seurat)
library(harmony)
library(AUCell)
library(ggpubr)
library(patchwork)
library(tidyverse)

setwd("D:/scRNA_project")

data_dir <- "D:/scRNA_project/GSE125449"

# --- Load Set1 and Set2 ---
set1_matrix <- ReadMtx(
  mtx      = file.path(data_dir, "GSE125449_Set1_matrix.mtx"),
  cells    = file.path(data_dir, "GSE125449_Set1_barcodes.tsv"),
  features = file.path(data_dir, "GSE125449_Set1_genes.tsv")
)
set2_matrix <- ReadMtx(
  mtx      = file.path(data_dir, "GSE125449_Set2_matrix.mtx"),
  cells    = file.path(data_dir, "GSE125449_Set2_barcodes.tsv"),
  features = file.path(data_dir, "GSE125449_Set2_genes.tsv")
)

seu_s1 <- CreateSeuratObject(set1_matrix, project = "Set1",
                             min.cells = 3, min.features = 200)
seu_s2 <- CreateSeuratObject(set2_matrix, project = "Set2",
                             min.cells = 3, min.features = 200)
seu_s1$batch <- "Set1"
seu_s2$batch <- "Set2"

seu_gse <- merge(seu_s1, y = seu_s2,
                 add.cell.ids = c("S1","S2"),
                 project = "GSE125449")

# --- QC ---
seu_gse[["percent.mt"]] <- PercentageFeatureSet(seu_gse, pattern = "^MT-")
seu_gse <- subset(seu_gse,
                  subset = nFeature_RNA > 200 &
                    nFeature_RNA < 6000 &
                    percent.mt < 20)
cat("Cells after QC:", ncol(seu_gse), "\n")

# --- Normalization, PCA, Harmony ---
seu_gse <- seu_gse %>%
  NormalizeData() %>%
  FindVariableFeatures(nfeatures = 3000) %>%
  ScaleData() %>%
  RunPCA(npcs = 50)

seu_gse <- RunHarmony(seu_gse,
                      group.by.vars = "batch",
                      dims.use = 1:30)
seu_gse <- seu_gse %>%
  RunUMAP(reduction = "harmony", dims = 1:30) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30) %>%
  FindClusters(resolution = 0.5)

# --- Extract hepatocytes (clusters 3, 5, 7 based on ALB/APOA1) ---
hepatocytes_gse <- subset(seu_gse, idents = c(3, 5, 7))
hepatocytes_gse <- JoinLayers(hepatocytes_gse)
cat("Hepatocyte cells:", ncol(hepatocytes_gse), "\n")

# --- AUCell glycolysis scoring ---
glycolysis_genes <- c(
  "HK1","HK2","GPI","PFKL","PFKP","PFKM",
  "ALDOA","ALDOB","ALDOC","TPI1","GAPDH",
  "PGK1","PGAM1","ENO1","ENO2","PKM","LDHA",
  "LDHB","SLC2A1","SLC2A3","PFKFB3","GCK"
)

expr_matrix <- GetAssayData(hepatocytes_gse, layer = "counts")
cells_rankings <- AUCell_buildRankings(expr_matrix,
                                       nCores = 1,
                                       plotStats = FALSE)
geneSets <- list(glycolysis = glycolysis_genes)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings,
                            aucMaxRank = ceiling(0.05 * nrow(cells_rankings)))

auc_scores <- as.numeric(getAUC(cells_AUC)["glycolysis", ])
hepatocytes_gse$glycolysis_score <- auc_scores

median_score <- median(auc_scores)
hepatocytes_gse$glyco_group <- ifelse(
  auc_scores >= median_score, "GlycoHigh", "GlycoLow"
)

# --- ENO1-glycolysis correlation ---
eno1_expr <- as.numeric(
  GetAssayData(hepatocytes_gse, layer = "data")["ENO1", ]
)
cor_result <- cor.test(eno1_expr, auc_scores, method = "pearson")
cat("ENO1-Glycolysis R =", round(cor_result$estimate, 3),
    ", p =", cor_result$p.value, "\n")

# --- SPP1 and MIF differential expression ---
spp1_expr <- as.numeric(
  GetAssayData(hepatocytes_gse, layer = "data")["SPP1", ]
)
mif_expr <- as.numeric(
  GetAssayData(hepatocytes_gse, layer = "data")["MIF", ]
)

p_spp1 <- ggboxplot(
  data.frame(SPP1 = spp1_expr,
             group = hepatocytes_gse$glyco_group),
  x = "group", y = "SPP1", fill = "group",
  palette = c("GlycoHigh" = "#E64B35",
              "GlycoLow"  = "#4DBBD5"),
  add = "jitter", xlab = "",
  ylab = "SPP1 Expression",
  title = "SPP1: GlycoHigh vs GlycoLow (GSE125449)"
) + stat_compare_means(method = "wilcox.test", label = "p.format")

p_mif <- ggboxplot(
  data.frame(MIF = mif_expr,
             group = hepatocytes_gse$glyco_group),
  x = "group", y = "MIF", fill = "group",
  palette = c("GlycoHigh" = "#E64B35",
              "GlycoLow"  = "#4DBBD5"),
  add = "jitter", xlab = "",
  ylab = "MIF Expression",
  title = "MIF: GlycoHigh vs GlycoLow (GSE125449)"
) + stat_compare_means(method = "wilcox.test", label = "p.format")

ggsave("FigS_GSE125449_SPP1_MIF_validation.pdf",
       p_spp1 + p_mif, width = 10, height = 5)

cat("GSE125449 validation complete.\n")
cat("Saved: FigS_GSE125449_SPP1_MIF_validation.pdf\n")