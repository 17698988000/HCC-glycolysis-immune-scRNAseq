# =============================================================================
# 05_GSE125449_validation.R
#
# Current purpose:
#   Independent GSE125449 source/QC-first validation of the glycolysis-hepatocyte
#   program without legacy local paths or obsolete cluster/label assumptions.
#
# This replacement removes:
#   - setwd("D:/scRNA_project")
#   - hard-coded D:/scRNA_project/GSE125449 paths
#   - fixed cluster IDs such as clusters 3, 5, and 7 as hepatocytes
#   - "Hepatocytes" plural as a required label
#
# Analysis performed:
#   1. Load GSE125449 from either a prebuilt Seurat RDS or Set1/Set2 matrix files.
#   2. Normalize and QC the object.
#   3. Identify hepatocyte-like cells using available annotation if present;
#      otherwise use lineage marker scoring.
#   4. Score glycolysis using AUCell and the locked glycolysis gene set.
#   5. Split hepatocyte-like cells into GlycoHigh / GlycoLow by the median
#      AUCell score, using GlycoHigh = score > median and GlycoLow = remaining.
#   6. Test ENO1-glycolysis correlation and SPP1/MIF/PTGES expression difference.
#
# Inputs searched automatically:
#   - GSE125449_seurat.rds
#   - data/GSE125449/GSE125449_seurat.rds
#   - GSE125449/GSE125449_seurat.rds
#   - data/GSE125449/GSE125449_Set1_matrix.mtx[.gz]
#   - data/GSE125449/GSE125449_Set1_barcodes.tsv[.gz]
#   - data/GSE125449/GSE125449_Set1_genes.tsv[.gz] or features.tsv[.gz]
#   - same for Set2
#
# Outputs:
#   results/GSE125449_05_qc.csv
#   results/GSE125449_05_hepatocyte_cell_metadata.csv
#   results/GSE125449_05_glycolysis_gene_presence.csv
#   results/GSE125449_05_validation_summary.csv
#   results/GSE125449_05_gene_group_tests.csv
#   results/GSE125449_05_lineage_marker_scores.csv
#   results/GSE125449_05_hepatocyte_scored.rds
#   results/GSE125449_05_validation_plot_no_title.pdf
#   results/GSE125449_05_validation_plot_no_title.png
#   results/GSE125449_05_manifest.csv
#   results/GSE125449_05_sessionInfo.txt
#
# Notes:
#   This script is a validation/source-data script, not a locked manuscript
#   figure-restoration script.
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(1234)

min_features <- 200
max_features <- 6000
max_percent_mt <- 20

glycolysis_genes <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKM",
  "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH",
  "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA",
  "LDHB", "SLC2A1", "SLC2A3", "PFKFB3", "GCK"
)

genes_to_test <- c("ENO1", "LDHA", "SPP1", "MIF", "PTGES")

lineage_markers <- list(
  Hepatocyte = c("ALB", "APOA1", "APOA2", "APOB", "TTR", "FGA", "FGB", "HP", "ASGR1", "CYP3A4"),
  T_NK = c("CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "KLRD1"),
  Myeloid = c("LYZ", "LST1", "C1QA", "C1QB", "C1QC", "S100A8", "S100A9"),
  B_cell = c("MS4A1", "CD79A", "CD79B", "MZB1"),
  Endothelial = c("PECAM1", "VWF", "KDR", "ENG"),
  Fibroblast = c("COL1A1", "COL1A2", "DCN", "LUM")
)

annotation_candidates <- c(
  "cell_type", "CellType", "celltype", "annotation", "Annotation",
  "cell_type_major", "major_cell_type", "majorType", "cell.types"
)

hepatocyte_annotation_pattern <- paste(
  c("hepatocyte", "hepatocytes", "malignant", "tumor", "hcc", "liver"),
  collapse = "|"
)

# -----------------------------------------------------------------------------
# Utilities
# -----------------------------------------------------------------------------

write_csv <- function(x, filename) {
  path <- file.path(out_dir, filename)
  utils::write.csv(x, path, row.names = FALSE, quote = TRUE)
  message("Wrote: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(path)
}

write_text <- function(lines, filename) {
  path <- file.path(out_dir, filename)
  writeLines(lines, con = path, useBytes = TRUE)
  message("Wrote: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(path)
}

required_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unknown"
  x
}

first_existing <- function(paths) {
  paths[file.exists(paths)][1]
}

find_existing_file <- function(root_dirs, file_names) {
  candidates <- as.vector(outer(root_dirs, file_names, file.path))
  first_existing(candidates)
}

get_assay_matrix <- function(seu, assay_use = NULL, layer_use = "data") {
  if (is.null(assay_use)) {
    assay_use <- Seurat::DefaultAssay(seu)
    if ("RNA" %in% Seurat::Assays(seu)) {
      assay_use <- "RNA"
    }
  }

  mat <- tryCatch(
    {
      if (requireNamespace("SeuratObject", quietly = TRUE) &&
          "LayerData" %in% getNamespaceExports("SeuratObject")) {
        SeuratObject::LayerData(seu, assay = assay_use, layer = layer_use)
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )

  if (is.null(mat)) {
    mat <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay_use, layer = layer_use),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    slot_use <- ifelse(layer_use == "counts", "counts", "data")
    mat <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay_use, slot = slot_use),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    stop("Could not retrieve assay matrix: assay=", assay_use, ", layer=", layer_use)
  }

  mat
}

mean_expr_by_gene <- function(expr_data, genes) {
  genes_present <- intersect(genes, rownames(expr_data))
  if (length(genes_present) == 0) {
    return(rep(NA_real_, ncol(expr_data)))
  }
  Matrix::colMeans(expr_data[genes_present, , drop = FALSE])
}

read_mtx_flexible <- function(mtx_file, cell_file, feature_file, project_name) {
  message("Reading matrix for ", project_name)
  message("  matrix:   ", mtx_file)
  message("  barcodes: ", cell_file)
  message("  features: ", feature_file)

  obj <- tryCatch(
    Seurat::ReadMtx(
      mtx = mtx_file,
      cells = cell_file,
      features = feature_file,
      feature.column = 2,
      unique.features = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(obj)) {
    obj <- Seurat::ReadMtx(
      mtx = mtx_file,
      cells = cell_file,
      features = feature_file,
      feature.column = 1,
      unique.features = TRUE
    )
  }

  obj
}

load_gse125449 <- function() {
  rds_candidates <- c(
    "GSE125449_seurat.rds",
    file.path("data", "GSE125449", "GSE125449_seurat.rds"),
    file.path("GSE125449", "GSE125449_seurat.rds"),
    file.path("results", "GSE125449_seurat.rds")
  )

  rds_file <- first_existing(rds_candidates)
  if (!is.na(rds_file)) {
    message("Loading prebuilt Seurat object: ", rds_file)
    seu <- readRDS(rds_file)
    seu$GSE125449_source <- "prebuilt_rds"
    return(seu)
  }

  root_dirs <- c(
    file.path("data", "GSE125449"),
    "GSE125449",
    "."
  )

  load_set <- function(set_name) {
    prefix <- paste0("GSE125449_", set_name)

    mtx <- find_existing_file(
      root_dirs,
      c(
        paste0(prefix, "_matrix.mtx"),
        paste0(prefix, "_matrix.mtx.gz"),
        paste0(set_name, "_matrix.mtx"),
        paste0(set_name, "_matrix.mtx.gz"),
        file.path(set_name, "matrix.mtx"),
        file.path(set_name, "matrix.mtx.gz")
      )
    )

    cells <- find_existing_file(
      root_dirs,
      c(
        paste0(prefix, "_barcodes.tsv"),
        paste0(prefix, "_barcodes.tsv.gz"),
        paste0(set_name, "_barcodes.tsv"),
        paste0(set_name, "_barcodes.tsv.gz"),
        file.path(set_name, "barcodes.tsv"),
        file.path(set_name, "barcodes.tsv.gz")
      )
    )

    features <- find_existing_file(
      root_dirs,
      c(
        paste0(prefix, "_genes.tsv"),
        paste0(prefix, "_genes.tsv.gz"),
        paste0(prefix, "_features.tsv"),
        paste0(prefix, "_features.tsv.gz"),
        paste0(set_name, "_genes.tsv"),
        paste0(set_name, "_genes.tsv.gz"),
        paste0(set_name, "_features.tsv"),
        paste0(set_name, "_features.tsv.gz"),
        file.path(set_name, "genes.tsv"),
        file.path(set_name, "genes.tsv.gz"),
        file.path(set_name, "features.tsv"),
        file.path(set_name, "features.tsv.gz")
      )
    )

    if (any(is.na(c(mtx, cells, features)))) {
      return(NULL)
    }

    mat <- read_mtx_flexible(mtx, cells, features, set_name)
    seu <- Seurat::CreateSeuratObject(
      counts = mat,
      project = set_name,
      min.cells = 3,
      min.features = min_features
    )
    seu$batch <- set_name
    seu$GSE125449_source <- "matrix_files"
    seu
  }

  seu_s1 <- load_set("Set1")
  seu_s2 <- load_set("Set2")

  if (is.null(seu_s1) && is.null(seu_s2)) {
    stop(
      "Could not find GSE125449 input. Provide either GSE125449_seurat.rds ",
      "or Set1/Set2 matrix/barcode/gene files under data/GSE125449/.",
      call. = FALSE
    )
  }

  if (!is.null(seu_s1) && !is.null(seu_s2)) {
    seu <- merge(
      seu_s1,
      y = seu_s2,
      add.cell.ids = c("Set1", "Set2"),
      project = "GSE125449"
    )
  } else if (!is.null(seu_s1)) {
    seu <- seu_s1
  } else {
    seu <- seu_s2
  }

  seu
}

# -----------------------------------------------------------------------------
# Package checks
# -----------------------------------------------------------------------------

required_package("Seurat")
required_package("AUCell")
required_package("Matrix")
required_package("ggplot2")
required_package("patchwork")

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
})

# -----------------------------------------------------------------------------
# Load, normalize, and QC object
# -----------------------------------------------------------------------------

seu <- load_gse125449()

if (!"batch" %in% colnames(seu@meta.data)) {
  if ("orig.ident" %in% colnames(seu@meta.data)) {
    seu$batch <- safe_chr(seu$orig.ident)
  } else {
    seu$batch <- "GSE125449"
  }
}

if (!"percent.mt" %in% colnames(seu@meta.data)) {
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
}

cells_before_qc <- ncol(seu)

seu <- subset(
  seu,
  subset = nFeature_RNA > min_features &
    nFeature_RNA < max_features &
    percent.mt < max_percent_mt
)

cells_after_qc <- ncol(seu)

if (cells_after_qc == 0) {
  stop("No cells retained after QC. Check input object and QC thresholds.", call. = FALSE)
}

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, verbose = FALSE)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Harmony is optional. It is used only if installed and more than one batch exists.
use_harmony <- requireNamespace("harmony", quietly = TRUE) && length(unique(seu$batch)) > 1

if (use_harmony) {
  message("Running Harmony batch correction.")
  seu <- harmony::RunHarmony(
    object = seu,
    group.by.vars = "batch",
    dims.use = 1:30,
    verbose = FALSE
  )
  reduction_use <- "harmony"
} else {
  message("Harmony not used; using PCA reduction.")
  reduction_use <- "pca"
}

seu <- RunUMAP(seu, reduction = reduction_use, dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = reduction_use, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

expr_data <- get_assay_matrix(seu, layer_use = "data")
expr_counts <- get_assay_matrix(seu, layer_use = "counts")

# -----------------------------------------------------------------------------
# Hepatocyte-like cell identification
# -----------------------------------------------------------------------------

annotation_col <- annotation_candidates[annotation_candidates %in% colnames(seu@meta.data)][1]
annotation_used <- !is.na(annotation_col)

if (annotation_used) {
  annot_values <- safe_chr(seu@meta.data[[annotation_col]])
  hepatocyte_like <- grepl(hepatocyte_annotation_pattern, annot_values, ignore.case = TRUE)
  hepatocyte_selection_method <- paste0("annotation:", annotation_col)
} else {
  lineage_scores <- sapply(
    lineage_markers,
    function(genes) mean_expr_by_gene(expr_data, genes)
  )

  lineage_scores <- as.data.frame(lineage_scores, stringsAsFactors = FALSE)
  lineage_scores$cell <- rownames(seu@meta.data)

  score_mat <- as.matrix(lineage_scores[, names(lineage_markers), drop = FALSE])
  max_lineage <- colnames(score_mat)[max.col(score_mat, ties.method = "first")]

  hepatocyte_like <- max_lineage == "Hepatocyte" &
    lineage_scores$Hepatocyte > 0

  hepatocyte_selection_method <- "lineage_marker_score"

  write_csv(lineage_scores, "GSE125449_05_lineage_marker_scores.csv")
}

hepatocyte_n <- sum(hepatocyte_like)

if (hepatocyte_n < 50) {
  stop(
    "Fewer than 50 hepatocyte-like cells were detected (n = ",
    hepatocyte_n,
    "). Review annotation/lineage marker selection before proceeding.",
    call. = FALSE
  )
}

hepa <- subset(seu, cells = colnames(seu)[hepatocyte_like])
expr_hepa_data <- get_assay_matrix(hepa, layer_use = "data")
expr_hepa_counts <- get_assay_matrix(hepa, layer_use = "counts")

# -----------------------------------------------------------------------------
# AUCell glycolysis scoring
# -----------------------------------------------------------------------------

gly_genes_present <- intersect(glycolysis_genes, rownames(expr_hepa_counts))
gly_gene_presence <- data.frame(
  gene = glycolysis_genes,
  present = glycolysis_genes %in% rownames(expr_hepa_counts),
  stringsAsFactors = FALSE
)

write_csv(gly_gene_presence, "GSE125449_05_glycolysis_gene_presence.csv")

if (length(gly_genes_present) < 8) {
  stop(
    "Too few glycolysis genes present for AUCell scoring: ",
    length(gly_genes_present),
    " of ",
    length(glycolysis_genes),
    call. = FALSE
  )
}

cells_rankings <- AUCell::AUCell_buildRankings(
  expr_hepa_counts,
  nCores = 1,
  plotStats = FALSE,
  verbose = FALSE
)

cells_AUC <- AUCell::AUCell_calcAUC(
  geneSets = list(glycolysis = gly_genes_present),
  rankings = cells_rankings,
  aucMaxRank = ceiling(0.05 * nrow(cells_rankings))
)

auc_scores <- as.numeric(AUCell::getAUC(cells_AUC)["glycolysis", ])
names(auc_scores) <- colnames(expr_hepa_counts)

median_score <- stats::median(auc_scores, na.rm = TRUE)
glyco_group <- ifelse(auc_scores > median_score, "GlycoHigh", "GlycoLow")

hepa$Glycolysis_AUC_GSE125449 <- auc_scores[colnames(hepa)]
hepa$glyco_group_GSE125449 <- glyco_group[colnames(hepa)]

glyco_high_n <- sum(hepa$glyco_group_GSE125449 == "GlycoHigh")
glyco_low_n <- sum(hepa$glyco_group_GSE125449 == "GlycoLow")

# -----------------------------------------------------------------------------
# Gene-level validation statistics
# -----------------------------------------------------------------------------

gene_stat_one <- function(gene) {
  if (!gene %in% rownames(expr_hepa_data)) {
    return(data.frame(
      gene = gene,
      present = FALSE,
      mean_GlycoHigh = NA_real_,
      mean_GlycoLow = NA_real_,
      mean_diff_High_minus_Low = NA_real_,
      median_GlycoHigh = NA_real_,
      median_GlycoLow = NA_real_,
      median_diff_High_minus_Low = NA_real_,
      wilcox_p = NA_real_,
      spearman_rho_vs_glycolysis = NA_real_,
      spearman_p_vs_glycolysis = NA_real_,
      pearson_r_vs_glycolysis = NA_real_,
      pearson_p_vs_glycolysis = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  x <- as.numeric(expr_hepa_data[gene, colnames(hepa)])
  group <- hepa$glyco_group_GSE125449

  high <- x[group == "GlycoHigh"]
  low <- x[group == "GlycoLow"]

  wilcox_p <- tryCatch(
    stats::wilcox.test(high, low)$p.value,
    error = function(e) NA_real_
  )

  sp <- tryCatch(
    suppressWarnings(stats::cor.test(x, auc_scores[colnames(hepa)], method = "spearman")),
    error = function(e) NULL
  )

  pr <- tryCatch(
    suppressWarnings(stats::cor.test(x, auc_scores[colnames(hepa)], method = "pearson")),
    error = function(e) NULL
  )

  data.frame(
    gene = gene,
    present = TRUE,
    mean_GlycoHigh = mean(high, na.rm = TRUE),
    mean_GlycoLow = mean(low, na.rm = TRUE),
    mean_diff_High_minus_Low = mean(high, na.rm = TRUE) - mean(low, na.rm = TRUE),
    median_GlycoHigh = stats::median(high, na.rm = TRUE),
    median_GlycoLow = stats::median(low, na.rm = TRUE),
    median_diff_High_minus_Low = stats::median(high, na.rm = TRUE) - stats::median(low, na.rm = TRUE),
    wilcox_p = wilcox_p,
    spearman_rho_vs_glycolysis = ifelse(is.null(sp), NA_real_, unname(sp$estimate)),
    spearman_p_vs_glycolysis = ifelse(is.null(sp), NA_real_, sp$p.value),
    pearson_r_vs_glycolysis = ifelse(is.null(pr), NA_real_, unname(pr$estimate)),
    pearson_p_vs_glycolysis = ifelse(is.null(pr), NA_real_, pr$p.value),
    stringsAsFactors = FALSE
  )
}

gene_tests <- do.call(rbind, lapply(genes_to_test, gene_stat_one))
gene_tests$wilcox_FDR <- stats::p.adjust(gene_tests$wilcox_p, method = "BH")
gene_tests$spearman_FDR <- stats::p.adjust(gene_tests$spearman_p_vs_glycolysis, method = "BH")
gene_tests$pearson_FDR <- stats::p.adjust(gene_tests$pearson_p_vs_glycolysis, method = "BH")

write_csv(gene_tests, "GSE125449_05_gene_group_tests.csv")

eno1_row <- gene_tests[gene_tests$gene == "ENO1", , drop = FALSE]

summary_df <- data.frame(
  item = c(
    "dataset",
    "input_source",
    "cells_before_qc",
    "cells_after_qc",
    "hepatocyte_selection_method",
    "annotation_column",
    "hepatocyte_like_cells",
    "glycolysis_genes_present",
    "glycolysis_genes_total",
    "median_Glycolysis_AUC",
    "GlycoHigh_n",
    "GlycoLow_n",
    "ENO1_pearson_r_vs_glycolysis",
    "ENO1_pearson_p_vs_glycolysis",
    "ENO1_spearman_rho_vs_glycolysis",
    "ENO1_spearman_p_vs_glycolysis"
  ),
  value = c(
    "GSE125449",
    paste(unique(seu$GSE125449_source), collapse = ";"),
    as.character(cells_before_qc),
    as.character(cells_after_qc),
    hepatocyte_selection_method,
    ifelse(annotation_used, annotation_col, "None"),
    as.character(hepatocyte_n),
    as.character(length(gly_genes_present)),
    as.character(length(glycolysis_genes)),
    sprintf("%.8f", median_score),
    as.character(glyco_high_n),
    as.character(glyco_low_n),
    ifelse(nrow(eno1_row) == 1, sprintf("%.6f", eno1_row$pearson_r_vs_glycolysis), NA),
    ifelse(nrow(eno1_row) == 1, format(eno1_row$pearson_p_vs_glycolysis, scientific = TRUE), NA),
    ifelse(nrow(eno1_row) == 1, sprintf("%.6f", eno1_row$spearman_rho_vs_glycolysis), NA),
    ifelse(nrow(eno1_row) == 1, format(eno1_row$spearman_p_vs_glycolysis, scientific = TRUE), NA)
  ),
  stringsAsFactors = FALSE
)

write_csv(summary_df, "GSE125449_05_validation_summary.csv")

qc_df <- data.frame(
  check_id = c(
    "no_legacy_setwd",
    "input_loaded",
    "qc_retained_cells",
    "hepatocyte_like_cells_detected",
    "glycolysis_gene_presence_minimum",
    "glycohigh_definition",
    "glycolow_definition",
    "glyco_groups_nonempty",
    "eno1_present",
    "spp1_present",
    "mif_present"
  ),
  expected = c(
    "No D:/scRNA_project setwd required",
    "GSE125449 object or matrix files",
    "> 0",
    ">= 50",
    ">= 8 glycolysis genes",
    "Glycolysis_AUC > median",
    "remaining hepatocyte-like cells",
    "Both GlycoHigh and GlycoLow present",
    "ENO1 present",
    "SPP1 present or recorded absent",
    "MIF present or recorded absent"
  ),
  observed = c(
    "No legacy setwd used",
    paste(unique(seu$GSE125449_source), collapse = ";"),
    as.character(cells_after_qc),
    as.character(hepatocyte_n),
    as.character(length(gly_genes_present)),
    "Glycolysis_AUC > median",
    "remaining hepatocyte-like cells",
    paste0("GlycoHigh=", glyco_high_n, "; GlycoLow=", glyco_low_n),
    as.character("ENO1" %in% rownames(expr_hepa_data)),
    as.character("SPP1" %in% rownames(expr_hepa_data)),
    as.character("MIF" %in% rownames(expr_hepa_data))
  ),
  pass = c(
    TRUE,
    TRUE,
    cells_after_qc > 0,
    hepatocyte_n >= 50,
    length(gly_genes_present) >= 8,
    TRUE,
    TRUE,
    glyco_high_n > 0 && glyco_low_n > 0,
    "ENO1" %in% rownames(expr_hepa_data),
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv(qc_df, "GSE125449_05_qc.csv")

if (!all(qc_df$pass)) {
  stop("Mandatory GSE125449 validation QC failed. See results/GSE125449_05_qc.csv.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Cell metadata output and plots
# -----------------------------------------------------------------------------

cell_meta <- data.frame(
  cell = colnames(hepa),
  batch = safe_chr(hepa$batch),
  seurat_cluster = safe_chr(hepa$seurat_clusters),
  hepatocyte_selection_method = hepatocyte_selection_method,
  Glycolysis_AUC_GSE125449 = hepa$Glycolysis_AUC_GSE125449,
  glyco_group_GSE125449 = hepa$glyco_group_GSE125449,
  stringsAsFactors = FALSE
)

if (annotation_used) {
  cell_meta$annotation_column <- annotation_col
  cell_meta$annotation_value <- safe_chr(hepa@meta.data[[annotation_col]])
}

for (gene in genes_to_test) {
  cell_meta[[gene]] <- if (gene %in% rownames(expr_hepa_data)) {
    as.numeric(expr_hepa_data[gene, colnames(hepa)])
  } else {
    NA_real_
  }
}

write_csv(cell_meta, "GSE125449_05_hepatocyte_cell_metadata.csv")

saveRDS(hepa, file.path(out_dir, "GSE125449_05_hepatocyte_scored.rds"))

plot_df <- cell_meta

p1 <- ggplot(plot_df, aes(x = Glycolysis_AUC_GSE125449, y = ENO1)) +
  geom_point(size = 0.25, alpha = 0.35) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  theme_classic(base_size = 10) +
  labs(
    x = "Glycolysis AUCell score",
    y = "ENO1 expression"
  )

plot_gene_box <- function(gene) {
  ggplot(plot_df, aes(x = glyco_group_GSE125449, y = .data[[gene]])) +
    geom_boxplot(outlier.size = 0.25, linewidth = 0.25) +
    theme_classic(base_size = 10) +
    labs(
      x = "",
      y = paste0(gene, " expression")
    )
}

p2 <- plot_gene_box("SPP1")
p3 <- plot_gene_box("MIF")
p4 <- plot_gene_box("PTGES")

combined_plot <- (p1 | p2) / (p3 | p4)

ggplot2::ggsave(
  filename = file.path(out_dir, "GSE125449_05_validation_plot_no_title.pdf"),
  plot = combined_plot,
  width = 8.5,
  height = 7
)

ggplot2::ggsave(
  filename = file.path(out_dir, "GSE125449_05_validation_plot_no_title.png"),
  plot = combined_plot,
  width = 8.5,
  height = 7,
  dpi = 300
)

# -----------------------------------------------------------------------------
# Manifest and session info
# -----------------------------------------------------------------------------

manifest <- data.frame(
  output_file = c(
    "GSE125449_05_qc.csv",
    "GSE125449_05_hepatocyte_cell_metadata.csv",
    "GSE125449_05_glycolysis_gene_presence.csv",
    "GSE125449_05_validation_summary.csv",
    "GSE125449_05_gene_group_tests.csv",
    "GSE125449_05_lineage_marker_scores.csv",
    "GSE125449_05_hepatocyte_scored.rds",
    "GSE125449_05_validation_plot_no_title.pdf",
    "GSE125449_05_validation_plot_no_title.png",
    "GSE125449_05_manifest.csv",
    "GSE125449_05_sessionInfo.txt"
  ),
  purpose = c(
    "Mandatory validation QC checks",
    "Cell-level hepatocyte-like metadata and expression source table",
    "Glycolysis gene presence table",
    "Dataset-level validation summary",
    "Gene-level group-difference and correlation statistics",
    "Lineage marker scores used when no annotation column is available",
    "Scored hepatocyte-like Seurat object",
    "No-title validation plot",
    "No-title validation plot",
    "Output manifest",
    "R session information"
  ),
  direct_manuscript_figure = c(
    "No", "No", "No", "No", "No", "No", "No",
    "No; validation plot only", "No; validation plot only", "No", "No"
  ),
  stringsAsFactors = FALSE
)

write_csv(manifest, "GSE125449_05_manifest.csv")
write_text(capture.output(sessionInfo()), "GSE125449_05_sessionInfo.txt")

message("")
message("GSE125449 validation complete.")
message("Cells before QC: ", cells_before_qc)
message("Cells after QC: ", cells_after_qc)
message("Hepatocyte-like cells: ", hepatocyte_n)
message("GlycoHigh: ", glyco_high_n)
message("GlycoLow: ", glyco_low_n)
message("ENO1 Pearson r vs glycolysis: ", eno1_row$pearson_r_vs_glycolysis)
message("Outputs written under: ", normalizePath(out_dir, winslash = "/", mustWork = FALSE))
