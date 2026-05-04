# ============================================================
# 06_inferCNV_malignant.R
# Purpose: Identify malignant hepatocytes using inferCNV and
#          transcriptomic malignancy marker scoring
# Corresponds to: Methods Section 2.5
# ============================================================

library(Seurat)
library(infercnv)
library(ggplot2)

# ── Load data ──────────────────────────────────────────────
seu <- readRDS("D:/scRNA_project/seurat_final.rds")

# ── 1. Transcriptomic malignancy marker scoring ────────────
malignancy_markers <- list(
  Malignancy = c("AFP", "GPC3", "EPCAM", "SALL4", "CDH17", "KRT19")
)

seu <- AddModuleScore(
  seu,
  features = malignancy_markers,
  name = "Malignancy_Score",
  ctrl = 100
)

# Classify hepatocytes above median as putatively malignant
hep_idx <- which(seu$cell_type == "Hepatocyte")
median_score <- median(seu$Malignancy_Score1[hep_idx])
seu$malignant_status <- "Non-Hepatocyte"
seu$malignant_status[hep_idx] <- ifelse(
  seu$Malignancy_Score1[hep_idx] >= median_score,
  "Malignant", "Non-malignant"
)

cat("Malignant hepatocytes:", sum(seu$malignant_status == "Malignant"), "\n")
cat("Non-malignant hepatocytes:", sum(seu$malignant_status == "Non-malignant"), "\n")

# ── 2. inferCNV ────────────────────────────────────────────
# Prepare count matrix (hepatocytes + reference T cells + endothelial cells)
ref_cells <- c("T cell", "Endothelial")
cells_use <- seu$cell_type %in% c("Hepatocyte", ref_cells)
seu_sub <- subset(seu, cells = colnames(seu)[cells_use])

counts_mat <- GetAssayData(seu_sub, assay = "RNA", layer = "counts")

# Cell annotations
cell_annots <- data.frame(
  cell_type = seu_sub$cell_type,
  row.names = colnames(seu_sub)
)
annot_file <- "D:/scRNA_project/inferCNV_output/cell_annotations.txt"
dir.create("D:/scRNA_project/inferCNV_output/", showWarnings = FALSE)
write.table(cell_annots, annot_file, sep = "\t", quote = FALSE, col.names = FALSE)

# Gene order file (hg38)
# Note: gene_order_file should be downloaded from inferCNV repository
# https://data.broadinstitute.org/Trinity/CTAT/cnv/hg38_gencode_v27.txt
gene_order_file <- "D:/scRNA_project/hg38_gencode_v27.txt"

infercnv_obj <- CreateInfercnvObject(
  raw_counts_matrix = counts_mat,
  annotations_file  = annot_file,
  delim             = "\t",
  gene_order_file   = gene_order_file,
  ref_group_names   = ref_cells
)

infercnv_obj <- infercnv::run(
  infercnv_obj,
  cutoff            = 0.1,
  out_dir           = "D:/scRNA_project/inferCNV_output/",
  cluster_by_groups = TRUE,
  denoise           = TRUE,
  HMM               = TRUE,
  HMM_type          = "i6",
  num_threads       = 4
)

cat("inferCNV analysis complete. Results saved to D:/scRNA_project/inferCNV_output/\n")
