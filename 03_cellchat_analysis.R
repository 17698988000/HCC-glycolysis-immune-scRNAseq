# ============================================================
# Script 03: Cell-Cell Communication Analysis (CellChat)
# ============================================================

library(Seurat)
library(CellChat)
library(tidyverse)

setwd("D:/scRNA_project")

seu <- readRDS("seurat_final.rds")

# --- Round 1: Overall communication (6 cell types) ---
data_input <- GetAssayData(seu, layer = "data")
meta <- data.frame(cell_type = seu$cell_type,
                   row.names = colnames(seu))

cellchat <- createCellChat(object = data_input,
                           meta = meta,
                           group.by = "cell_type")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "hcc_cellchat.rds")

# --- Round 2: Glycolysis-stratified (7 groups) ---
# Add glyco_group to metadata
seu$cell_type_glyco <- seu$cell_type
hepa_cells <- colnames(seu)[seu$cell_type == "Hepatocytes" & 
                              seu$site == "Tumor"]

tumor_hepa <- readRDS("tumor_hepatocytes.rds")
seu$cell_type_glyco[hepa_cells] <- tumor_hepa$glyco_group[hepa_cells]

meta_glyco <- data.frame(
  cell_type = seu$cell_type_glyco,
  row.names = colnames(seu)
)

cellchat2 <- createCellChat(object = GetAssayData(seu, layer = "data"),
                            meta = meta_glyco,
                            group.by = "cell_type")
cellchat2@DB <- CellChatDB.human

cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2, type = "triMean",
                               nboot = 1000)  # permutation test
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)

saveRDS(cellchat2, "hcc_cellchat2_glyco.rds")
cat("CellChat analysis complete.\n")