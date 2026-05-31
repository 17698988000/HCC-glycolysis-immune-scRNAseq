# ============================================================
# Script 01: Quality Control, Clustering, and Cell Annotation
# Dataset: GSE149614 (8 HCC patients, 54,631 cells)
# ============================================================

library(Seurat)
library(harmony)
library(tidyverse)

# Run from the repository root. Set PROJECT_DIR before running only when needed.
project_dir <- Sys.getenv("PROJECT_DIR", unset = ".")
if (!identical(project_dir, ".")) setwd(project_dir)

# --- Load raw count matrix ---
# Data downloaded from GEO: GSE149614
# Subset to 8 patients with matched PT and NTL samples (HCC03-HCC10)
counts <- Read10X("GSE149614_HCC.scRNAseq.S71915.count.txt.gz")
metadata <- read.table("GSE149614_HCC.metadata.txt.gz", 
                       header = TRUE, sep = "\t")

# --- Create Seurat object ---
seu <- CreateSeuratObject(counts = counts,
                          meta.data = metadata,
                          min.cells = 3,
                          min.features = 200)

# --- Quality control ---
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- subset(seu, subset = 
                nFeature_RNA > 200 &
                nFeature_RNA < 6000 &
                nCount_RNA > 500 &
                percent.mt < 20)

cat("Cells after QC:", ncol(seu), "\n")

# --- Subset to 8 patients with matched PT+NTL ---
patients_keep <- c("HCC03","HCC04","HCC05","HCC06",
                   "HCC07","HCC08","HCC09","HCC10")
seu <- subset(seu, subset = patient %in% patients_keep)
cat("Cells after patient subsetting:", ncol(seu), "\n")

# --- Normalization and HVG ---
seu <- NormalizeData(seu, normalization.method = "LogNormalize",
                     scale.factor = 10000)
seu <- FindVariableFeatures(seu, nfeatures = 3000)

# --- PCA ---
seu <- ScaleData(seu)
seu <- RunPCA(seu, npcs = 50)

# --- Harmony batch correction (by patient) ---
seu <- RunHarmony(seu, group.by.vars = "patient", dims.use = 1:30)

# --- Clustering and UMAP ---
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30)

# --- Cell type annotation ---
# Marker genes for each cell type
markers <- list(
  B_cells       = c("CD79A", "MS4A1"),
  Endothelial   = c("PECAM1", "VWF"),
  Fibroblasts   = c("COL1A1", "DCN"),
  Hepatocytes   = c("ALB", "APOA1"),
  Myeloid       = c("CD68", "CD14"),
  TNK_cells     = c("CD3D", "NKG7"),
  CD8T          = c("CD8A"),
  LSEC          = c("CLEC4M"),
  pDC           = c("CLEC4C", "LILRA4")
)

# Assign cell types based on marker expression
# (Manual annotation based on DotPlot and VlnPlot inspection)
cell_type_map <- c(
  "0"  = "Myeloid",
  "1"  = "T_NK",
  "2"  = "Hepatocytes",
  "3"  = "T_NK",
  "4"  = "Myeloid",
  "5"  = "Hepatocytes",
  "6"  = "B_cells",
  "7"  = "Endothelial",
  "8"  = "Fibroblasts",
  "9"  = "CD8T",
  "10" = "LSEC",
  "11" = "pDC"
)
seu$cell_type <- cell_type_map[as.character(seu$seurat_clusters)]

# --- Save ---
saveRDS(seu, "seurat_8patients.rds")
cat("Saved: seurat_8patients.rds\n")