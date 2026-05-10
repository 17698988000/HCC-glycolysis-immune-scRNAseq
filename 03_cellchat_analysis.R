# ============================================================
# Script 03: Cell-Cell Communication Analysis (CellChat)
# Section 2.6
# Round 1: 6 original cell types (overall landscape)
# Round 2: 7 groups (GlycoHigh/GlycoLow hepatocytes + 5 immune)
# Round 3: Subtype-resolved (T/NK and myeloid subclustered)
# ============================================================

library(Seurat)
library(CellChat)
library(tidyverse)

setwd("D:/scRNA_project")

seu        <- readRDS("seurat_final.rds")
tumor_hepa <- readRDS("tumor_hepatocytes.rds")

# ============================================================
# Round 1: Overall communication using 6 original cell types
# ============================================================
cat("=== Round 1: Overall CellChat (6 cell types) ===\n")

data_input <- GetAssayData(seu, layer = "data")
meta_r1    <- data.frame(cell_type = seu$cell_type,
                         row.names = colnames(seu))

cellchat <- createCellChat(object = data_input,
                           meta   = meta_r1,
                           group.by = "cell_type")
cellchat@DB <- CellChatDB.human

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

saveRDS(cellchat, "hcc_cellchat_round1.rds")
cat("Round 1 complete. Saved: hcc_cellchat_round1.rds\n")

# ============================================================
# Round 2: Glycolysis-stratified (7 groups)
# TumorGlycoHigh / TumorGlycoLow + 5 other cell types
# ============================================================
cat("\n=== Round 2: Glycolysis-stratified CellChat (7 groups) ===\n")

# Build cell_type_glyco: replace tumor hepatocytes with GlycoHigh/GlycoLow
seu$cell_type_glyco <- seu$cell_type

hepa_cells <- colnames(seu)[seu$cell_type == "Hepatocytes" &
                               seu$site    == "Tumor"]
cat("Tumor hepatocytes to recode:", length(hepa_cells), "\n")

# Align barcodes (tumor_hepa may be a subset of seu)
shared_cells <- intersect(hepa_cells, colnames(tumor_hepa))
seu$cell_type_glyco[shared_cells] <- tumor_hepa$glyco_group[shared_cells]

cat("GlycoHigh cells:", sum(seu$cell_type_glyco == "TumorGlycoHigh"), "\n")
cat("GlycoLow  cells:", sum(seu$cell_type_glyco == "TumorGlycoLow"),  "\n")

# ── BUG FIX: group.by must reference the glyco-stratified column ──
meta_r2 <- data.frame(cell_type_glyco = seu$cell_type_glyco,
                      row.names       = colnames(seu))

cellchat2 <- createCellChat(object   = GetAssayData(seu, layer = "data"),
                            meta     = meta_r2,
                            group.by = "cell_type_glyco")   # ← FIXED (was "cell_type")
cellchat2@DB <- CellChatDB.human

cellchat2 <- subsetData(cellchat2)
cellchat2 <- identifyOverExpressedGenes(cellchat2)
cellchat2 <- identifyOverExpressedInteractions(cellchat2)
cellchat2 <- computeCommunProb(cellchat2, type  = "triMean",
                                           nboot = 1000)   # permutation test
cellchat2 <- filterCommunication(cellchat2, min.cells = 10)
cellchat2 <- computeCommunProbPathway(cellchat2)
cellchat2 <- aggregateNet(cellchat2)

saveRDS(cellchat2, "hcc_cellchat_round2_glyco.rds")
cat("Round 2 complete. Saved: hcc_cellchat_round2_glyco.rds\n")

# ============================================================
# Round 3: Subtype-resolved CellChat
# T/NK and myeloid cells subclustered (Section 2.6, Supplementary Figure S5)
# ============================================================
cat("\n=== Round 3: Subtype-resolved CellChat ===\n")

# ── Step 1: Subcluster T/NK cells ──
tnk_cells <- subset(seu, cell_type %in% c("T_NK", "CD8T"))
tnk_cells <- FindNeighbors(tnk_cells, reduction = "harmony", dims = 1:15)
tnk_cells <- FindClusters(tnk_cells, resolution = 0.4)

# Annotate T/NK subtypes based on marker expression
# Exhausted T:      TIGIT, HAVCR2, TOX
# Effector memory:  NKG7, CCR7
# Activated T:      PDCD1, GZMB
# NK cells:         KLRD1, GNLY
# Naive T:          TCF7, SELL
tnk_markers <- list(
  TNK_Tex = c("TIGIT", "HAVCR2", "TOX"),
  TNK_Tem = c("NKG7",  "CCR7"),
  TNK_Tact= c("PDCD1", "GZMB"),
  TNK_NK  = c("KLRD1", "GNLY"),
  TNK_Tn  = c("TCF7",  "SELL")
)

# Score each cluster and assign subtype by highest mean score
tnk_expr    <- GetAssayData(tnk_cells, layer = "data")
cluster_ids <- levels(tnk_cells$seurat_clusters)

tnk_subtype_map <- sapply(cluster_ids, function(cl) {
  cl_cells <- colnames(tnk_cells)[tnk_cells$seurat_clusters == cl]
  scores   <- sapply(tnk_markers, function(genes) {
    g <- intersect(genes, rownames(tnk_expr))
    if (length(g) == 0) return(0)
    mean(colMeans(tnk_expr[g, cl_cells, drop = FALSE]))
  })
  names(which.max(scores))
})

tnk_cells$cell_subtype <- tnk_subtype_map[as.character(tnk_cells$seurat_clusters)]
cat("T/NK subtypes:\n"); print(table(tnk_cells$cell_subtype))

# ── Step 2: Subcluster myeloid cells ──
mye_cells <- subset(seu, cell_type == "Myeloid")
mye_cells <- FindNeighbors(mye_cells, reduction = "harmony", dims = 1:15)
mye_cells <- FindClusters(mye_cells, resolution = 0.4)

# Kupffer cells:       C1QA, APOE, FOLR2
# Tumor-assoc. macro:  CD163, CCL18
# MDSC:                S100A9, LYZ, VCAN
# Monocyte-derived:    HLA-DRA, LYZ
mye_markers <- list(
  Mye_KC   = c("C1QA",  "APOE",  "FOLR2"),
  Mye_TAM  = c("CD163", "CCL18"),
  Mye_MDSC = c("S100A9","LYZ",   "VCAN"),
  Mye_Mono = c("HLA-DRA","LYZ")
)

mye_expr        <- GetAssayData(mye_cells, layer = "data")
mye_cluster_ids <- levels(mye_cells$seurat_clusters)

mye_subtype_map <- sapply(mye_cluster_ids, function(cl) {
  cl_cells <- colnames(mye_cells)[mye_cells$seurat_clusters == cl]
  scores   <- sapply(mye_markers, function(genes) {
    g <- intersect(genes, rownames(mye_expr))
    if (length(g) == 0) return(0)
    mean(colMeans(mye_expr[g, cl_cells, drop = FALSE]))
  })
  names(which.max(scores))
})

mye_cells$cell_subtype <- mye_subtype_map[as.character(mye_cells$seurat_clusters)]
cat("Myeloid subtypes:\n"); print(table(mye_cells$cell_subtype))

# ── Step 3: Build refined annotation for Round 3 ──
seu$cell_subtype <- seu$cell_type_glyco   # start from Round 2 annotation

# Overwrite T/NK and myeloid with subtype labels
tnk_shared <- intersect(colnames(tnk_cells), colnames(seu))
mye_shared <- intersect(colnames(mye_cells), colnames(seu))
seu$cell_subtype[tnk_shared] <- tnk_cells$cell_subtype[tnk_shared]
seu$cell_subtype[mye_shared] <- mye_cells$cell_subtype[mye_shared]

cat("\nRound 3 cell type distribution:\n")
print(table(seu$cell_subtype))

# ── Step 4: Run CellChat Round 3 ──
meta_r3 <- data.frame(cell_subtype = seu$cell_subtype,
                      row.names    = colnames(seu))

cellchat3 <- createCellChat(object   = GetAssayData(seu, layer = "data"),
                            meta     = meta_r3,
                            group.by = "cell_subtype")
cellchat3@DB <- CellChatDB.human

cellchat3 <- subsetData(cellchat3)
cellchat3 <- identifyOverExpressedGenes(cellchat3)
cellchat3 <- identifyOverExpressedInteractions(cellchat3)
cellchat3 <- computeCommunProb(cellchat3, type  = "triMean",
                                           nboot = 1000)
cellchat3 <- filterCommunication(cellchat3, min.cells = 10)
cellchat3 <- computeCommunProbPathway(cellchat3)
cellchat3 <- aggregateNet(cellchat3)

saveRDS(cellchat3, "hcc_cellchat_round3_subtypes.rds")
cat("Round 3 complete. Saved: hcc_cellchat_round3_subtypes.rds\n")

cat("\n=== All CellChat analyses complete ===\n")
cat("Outputs:\n")
cat("  hcc_cellchat_round1.rds        - Overall 6-cell-type communication\n")
cat("  hcc_cellchat_round2_glyco.rds  - GlycoHigh/GlycoLow stratified\n")
cat("  hcc_cellchat_round3_subtypes.rds - Subtype-resolved\n")
