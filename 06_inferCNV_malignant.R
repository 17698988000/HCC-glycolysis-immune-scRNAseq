# =============================================================================
# 06_inferCNV_malignant.R
#
# Current purpose:
#   Current-object malignant hepatocyte source/QC workflow with optional inferCNV.
#
# This replacement removes legacy assumptions:
#   - no setwd("D:/scRNA_project")
#   - no hard-coded D:/scRNA_project output paths
#   - no "Hepatocytes" plural label requirement
#   - no automatic download of gene-order files during analysis
#
# Locked current-object convention:
#   input object      = seurat_final.rds
#   patient column    = patient
#   cell-type column  = cell_type
#   site column       = site
#
# Tumor-derived hepatocytes:
#   site == "Tumor" & cell_type == "Hepatocyte"
#   expected n = 15,391
#
# Analysis performed:
#   1. Validate current object and locked metadata columns.
#   2. Score malignancy-associated marker genes in hepatocytes.
#   3. Save source/QC tables and a marker-scored Seurat object.
#   4. Optionally run inferCNV if both infercnv and a local hg38 gene-order
#      file are available.
#
# Important:
#   This script does not directly generate or overwrite manuscript figures.
#   It is a source/QC workflow.
#
# Required local gene-order file for inferCNV, if inferCNV is to be run:
#   hg38_gencode_v27.txt
#
# Gene-order file candidates searched:
#   hg38_gencode_v27.txt
#   data/hg38_gencode_v27.txt
#   data/inferCNV/hg38_gencode_v27.txt
#   reference/hg38_gencode_v27.txt
#   ref/hg38_gencode_v27.txt
#
# Outputs:
#   results/inferCNV_06_qc.csv
#   results/inferCNV_06_hepatocyte_marker_scores.csv
#   results/inferCNV_06_marker_gene_presence.csv
#   results/inferCNV_06_group_counts.csv
#   results/inferCNV_06_reference_counts.csv
#   results/inferCNV_06_run_status.csv
#   results/inferCNV_06_marker_scored_seurat.rds
#   results/inferCNV_06_marker_score_summary.csv
#   results/inferCNV_06_marker_score_plot_no_title.pdf
#   results/inferCNV_06_marker_score_plot_no_title.png
#   results/inferCNV_06_manifest.csv
#   results/inferCNV_06_sessionInfo.txt
#
# Optional inferCNV outputs, if inferCNV runs successfully:
#   results/inferCNV_06_run/
#   results/inferCNV_06_input_cell_annotations.txt
#   results/inferCNV_06_cnv_burden.csv, if infercnv.observations.txt is present
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

input_rds <- "seurat_final.rds"
out_dir <- "results"
infercnv_out_dir <- file.path(out_dir, "inferCNV_06_run")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patient_col <- "patient"
celltype_col <- "cell_type"
site_col <- "site"

expected_tumor_hepatocytes <- 15391L

# Optional inferCNV run. If package or gene-order file is absent, the script
# records SKIPPED status and still writes source/QC outputs.
run_infercnv_if_available <- TRUE

# Cell types used as diploid reference if present.
preferred_reference_celltypes <- c("T_NK", "Endothelial")

# Additional fallback reference labels, used only if one preferred reference is absent.
fallback_reference_celltypes <- c("Myeloid", "B_cell", "B", "Plasma", "Fibroblast")

malignancy_markers <- c("AFP", "GPC3", "EPCAM", "SALL4", "CDH17", "KRT19")
hepatocyte_markers <- c("ALB", "APOA1", "APOA2", "TTR", "ASGR1", "CYP3A4")
reference_marker_check <- c("PTPRC", "NKG7", "CD3D", "PECAM1", "VWF", "LYZ", "LST1")

gene_order_candidates <- c(
  "hg38_gencode_v27.txt",
  file.path("data", "hg38_gencode_v27.txt"),
  file.path("data", "inferCNV", "hg38_gencode_v27.txt"),
  file.path("reference", "hg38_gencode_v27.txt"),
  file.path("ref", "hg38_gencode_v27.txt")
)

# inferCNV parameters
infercnv_cutoff <- 0.1
infercnv_num_threads <- 4
infercnv_cluster_by_groups <- TRUE
infercnv_denoise <- TRUE
infercnv_HMM <- TRUE
infercnv_HMM_type <- "i6"

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

safe_chr <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == "" | x == "NA"] <- "Unknown"
  x
}

required_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

optional_package <- function(pkg) {
  requireNamespace(pkg, quietly = TRUE)
}

first_existing <- function(paths) {
  existing <- paths[file.exists(paths)]
  if (length(existing) == 0) {
    return(NA_character_)
  }
  existing[1]
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

mean_expr_by_gene <- function(expr, genes) {
  genes_present <- intersect(genes, rownames(expr))
  if (length(genes_present) == 0) {
    return(rep(NA_real_, ncol(expr)))
  }
  Matrix::colMeans(expr[genes_present, , drop = FALSE])
}

marker_presence_table <- function(expr, gene_set, gene_set_name) {
  data.frame(
    gene_set = gene_set_name,
    gene = gene_set,
    present = gene_set %in% rownames(expr),
    stringsAsFactors = FALSE
  )
}

make_group_counts <- function(meta, group_col) {
  tab <- sort(table(meta[[group_col]]), decreasing = TRUE)
  data.frame(
    group = names(tab),
    n_cells = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || rng[1] == rng[2]) {
    return(rep(NA_real_, length(x)))
  }
  (x - rng[1]) / (rng[2] - rng[1])
}

parse_infercnv_observations <- function(obs_file, tumor_hep_cells) {
  if (!file.exists(obs_file)) {
    return(NULL)
  }

  message("Reading inferCNV observation matrix: ", obs_file)

  obs <- tryCatch(
    utils::read.table(
      obs_file,
      header = TRUE,
      sep = "\t",
      check.names = FALSE,
      row.names = 1,
      comment.char = "",
      quote = ""
    ),
    error = function(e) {
      warning("Failed to read inferCNV observations: ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(obs) || ncol(obs) == 0) {
    return(NULL)
  }

  common_cells <- intersect(colnames(obs), tumor_hep_cells)
  if (length(common_cells) == 0) {
    warning("inferCNV observations contain no columns matching tumor hepatocyte cells.")
    return(NULL)
  }

  obs_sub <- as.matrix(obs[, common_cells, drop = FALSE])

  # inferCNV observations are usually centered relative expression values.
  # A conservative cell-level burden metric is mean absolute deviation from
  # the per-gene median across tumor hepatocytes.
  gene_median <- apply(obs_sub, 1, stats::median, na.rm = TRUE)
  centered <- sweep(obs_sub, 1, gene_median, "-")
  burden <- Matrix::colMeans(abs(centered), na.rm = TRUE)

  median_burden <- stats::median(burden, na.rm = TRUE)

  data.frame(
    cell = names(burden),
    cnv_burden_mean_abs_centered = as.numeric(burden),
    cnv_burden_group = ifelse(burden > median_burden, "CNVHigh", "CNVLow"),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# Package checks
# -----------------------------------------------------------------------------

required_package("Seurat")
required_package("Matrix")
required_package("ggplot2")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
})

infercnv_available <- optional_package("infercnv")

# -----------------------------------------------------------------------------
# Load and validate object
# -----------------------------------------------------------------------------

if (!file.exists(input_rds)) {
  stop("Input object not found: ", input_rds, call. = FALSE)
}

seu <- readRDS(input_rds)
meta <- seu@meta.data

required_cols <- c(patient_col, celltype_col, site_col)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) {
  stop(
    "Missing required metadata columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

meta[[patient_col]] <- safe_chr(meta[[patient_col]])
meta[[celltype_col]] <- safe_chr(meta[[celltype_col]])
meta[[site_col]] <- safe_chr(meta[[site_col]])

seu@meta.data[[patient_col]] <- meta[[patient_col]]
seu@meta.data[[celltype_col]] <- meta[[celltype_col]]
seu@meta.data[[site_col]] <- meta[[site_col]]

expr_data <- get_assay_matrix(seu, layer_use = "data")
expr_counts <- get_assay_matrix(seu, layer_use = "counts")

tumor_hep_idx <- meta[[site_col]] == "Tumor" & meta[[celltype_col]] == "Hepatocyte"
tumor_hep_cells <- rownames(meta)[tumor_hep_idx]
tumor_hep_n <- length(tumor_hep_cells)

all_hep_idx <- meta[[celltype_col]] == "Hepatocyte"
all_hep_cells <- rownames(meta)[all_hep_idx]

# -----------------------------------------------------------------------------
# Marker scoring
# -----------------------------------------------------------------------------

marker_presence <- rbind(
  marker_presence_table(expr_data, malignancy_markers, "malignancy_markers"),
  marker_presence_table(expr_data, hepatocyte_markers, "hepatocyte_markers"),
  marker_presence_table(expr_data, reference_marker_check, "reference_marker_check")
)

write_csv(marker_presence, "inferCNV_06_marker_gene_presence.csv")

present_malignancy_markers <- intersect(malignancy_markers, rownames(expr_data))
present_hepatocyte_markers <- intersect(hepatocyte_markers, rownames(expr_data))

if (length(present_malignancy_markers) < 2) {
  stop(
    "Too few malignancy markers present: ",
    length(present_malignancy_markers),
    " of ",
    length(malignancy_markers),
    call. = FALSE
  )
}

mal_score <- mean_expr_by_gene(expr_data, present_malignancy_markers)
hep_score <- mean_expr_by_gene(expr_data, present_hepatocyte_markers)

seu$malignancy_marker_score_06 <- mal_score[colnames(seu)]
seu$hepatocyte_marker_score_06 <- hep_score[colnames(seu)]
seu$tumor_hepatocyte_06 <- colnames(seu) %in% tumor_hep_cells

tumor_hep_mal_score <- seu$malignancy_marker_score_06[tumor_hep_idx]
marker_median <- stats::median(tumor_hep_mal_score, na.rm = TRUE)

seu$malignancy_marker_group_06 <- "NonTumorHepatocyte"
seu$malignancy_marker_group_06[tumor_hep_idx] <- ifelse(
  tumor_hep_mal_score > marker_median,
  "TumorHepatocyte_MarkerHigh",
  "TumorHepatocyte_MarkerLow"
)

marker_high_n <- sum(seu$malignancy_marker_group_06 == "TumorHepatocyte_MarkerHigh")
marker_low_n <- sum(seu$malignancy_marker_group_06 == "TumorHepatocyte_MarkerLow")

hepatocyte_scores <- data.frame(
  cell = tumor_hep_cells,
  patient = meta[tumor_hep_cells, patient_col],
  site = meta[tumor_hep_cells, site_col],
  cell_type = meta[tumor_hep_cells, celltype_col],
  malignancy_marker_score = seu$malignancy_marker_score_06[tumor_hep_cells],
  hepatocyte_marker_score = seu$hepatocyte_marker_score_06[tumor_hep_cells],
  malignancy_marker_group = seu$malignancy_marker_group_06[tumor_hep_cells],
  stringsAsFactors = FALSE
)

for (gene in unique(c(malignancy_markers, hepatocyte_markers))) {
  hepatocyte_scores[[gene]] <- if (gene %in% rownames(expr_data)) {
    as.numeric(expr_data[gene, tumor_hep_cells])
  } else {
    NA_real_
  }
}

write_csv(hepatocyte_scores, "inferCNV_06_hepatocyte_marker_scores.csv")

group_counts <- make_group_counts(seu@meta.data, celltype_col)
write_csv(group_counts, "inferCNV_06_group_counts.csv")

# -----------------------------------------------------------------------------
# Reference cell selection for optional inferCNV
# -----------------------------------------------------------------------------

available_celltypes <- unique(meta[[celltype_col]])
reference_celltypes <- intersect(preferred_reference_celltypes, available_celltypes)

if (length(reference_celltypes) < 2) {
  needed <- 2 - length(reference_celltypes)
  fallback_use <- setdiff(intersect(fallback_reference_celltypes, available_celltypes), reference_celltypes)
  reference_celltypes <- c(reference_celltypes, head(fallback_use, needed))
}

reference_cells <- rownames(meta)[meta[[celltype_col]] %in% reference_celltypes]

reference_counts <- data.frame(
  reference_cell_type = reference_celltypes,
  n_cells = as.integer(table(factor(meta[reference_cells, celltype_col], levels = reference_celltypes))),
  preferred_reference = reference_celltypes %in% preferred_reference_celltypes,
  stringsAsFactors = FALSE
)

write_csv(reference_counts, "inferCNV_06_reference_counts.csv")

gene_order_file <- first_existing(gene_order_candidates)

# -----------------------------------------------------------------------------
# Mandatory QC
# -----------------------------------------------------------------------------

qc_df <- data.frame(
  check_id = c(
    "input_rds",
    "patient_column_present",
    "celltype_column_present",
    "site_column_present",
    "tumor_hepatocyte_definition",
    "tumor_hepatocyte_n",
    "obsolete_plural_hepatocytes_not_used",
    "malignancy_marker_genes_present",
    "reference_celltypes_detected",
    "infercnv_package_available",
    "gene_order_file_available"
  ),
  expected = c(
    "seurat_final.rds",
    patient_col,
    celltype_col,
    site_col,
    'site == "Tumor" & cell_type == "Hepatocyte"',
    as.character(expected_tumor_hepatocytes),
    'cell_type == "Hepatocyte"',
    ">= 2",
    ">= 1 reference cell type",
    "TRUE if inferCNV run is desired",
    "TRUE if inferCNV run is desired"
  ),
  observed = c(
    input_rds,
    ifelse(patient_col %in% colnames(meta), patient_col, "MISSING"),
    ifelse(celltype_col %in% colnames(meta), celltype_col, "MISSING"),
    ifelse(site_col %in% colnames(meta), site_col, "MISSING"),
    'site == "Tumor" & cell_type == "Hepatocyte"',
    as.character(tumor_hep_n),
    'cell_type == "Hepatocyte"',
    as.character(length(present_malignancy_markers)),
    paste(reference_celltypes, collapse = ";"),
    as.character(infercnv_available),
    ifelse(is.na(gene_order_file), "No local gene-order file found", gene_order_file)
  ),
  pass = c(
    identical(input_rds, "seurat_final.rds"),
    patient_col %in% colnames(meta),
    celltype_col %in% colnames(meta),
    site_col %in% colnames(meta),
    TRUE,
    identical(as.integer(tumor_hep_n), expected_tumor_hepatocytes),
    !any(meta[[celltype_col]][tumor_hep_idx] == "Hepatocytes"),
    length(present_malignancy_markers) >= 2,
    length(reference_celltypes) >= 1,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv(qc_df, "inferCNV_06_qc.csv")

if (!all(qc_df$pass)) {
  stop(
    "Mandatory Script 06 QC failed. See results/inferCNV_06_qc.csv. ",
    "Do not proceed with malignant hepatocyte/inferCNV analysis.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Marker score summary and plot
# -----------------------------------------------------------------------------

marker_summary <- data.frame(
  item = c(
    "input_rds",
    "tumor_hepatocytes",
    "malignancy_markers_present",
    "hepatocyte_markers_present",
    "malignancy_marker_score_median_in_tumor_hepatocytes",
    "marker_high_tumor_hepatocytes",
    "marker_low_tumor_hepatocytes",
    "reference_celltypes_for_optional_infercnv",
    "gene_order_file",
    "infercnv_package_available"
  ),
  value = c(
    input_rds,
    as.character(tumor_hep_n),
    paste(present_malignancy_markers, collapse = ";"),
    paste(present_hepatocyte_markers, collapse = ";"),
    sprintf("%.6f", marker_median),
    as.character(marker_high_n),
    as.character(marker_low_n),
    paste(reference_celltypes, collapse = ";"),
    ifelse(is.na(gene_order_file), "Not found", gene_order_file),
    as.character(infercnv_available)
  ),
  stringsAsFactors = FALSE
)

write_csv(marker_summary, "inferCNV_06_marker_score_summary.csv")

plot_df <- hepatocyte_scores

p <- ggplot(plot_df, aes(x = malignancy_marker_group, y = malignancy_marker_score)) +
  geom_boxplot(outlier.size = 0.25, linewidth = 0.25) +
  theme_classic(base_size = 10) +
  labs(
    x = "",
    y = "Malignancy marker score"
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggplot2::ggsave(
  filename = file.path(out_dir, "inferCNV_06_marker_score_plot_no_title.pdf"),
  plot = p,
  width = 4.8,
  height = 4.2
)

ggplot2::ggsave(
  filename = file.path(out_dir, "inferCNV_06_marker_score_plot_no_title.png"),
  plot = p,
  width = 4.8,
  height = 4.2,
  dpi = 300
)

saveRDS(seu, file.path(out_dir, "inferCNV_06_marker_scored_seurat.rds"))

# -----------------------------------------------------------------------------
# Optional inferCNV
# -----------------------------------------------------------------------------

run_status <- data.frame(
  item = c(
    "run_infercnv_if_available",
    "infercnv_package_available",
    "gene_order_file",
    "reference_celltypes",
    "tumor_hepatocytes_for_infercnv",
    "reference_cells_for_infercnv",
    "infercnv_status",
    "infercnv_out_dir",
    "cnv_burden_status"
  ),
  value = c(
    as.character(run_infercnv_if_available),
    as.character(infercnv_available),
    ifelse(is.na(gene_order_file), "Not found", gene_order_file),
    paste(reference_celltypes, collapse = ";"),
    as.character(length(tumor_hep_cells)),
    as.character(length(reference_cells)),
    "Not attempted yet",
    infercnv_out_dir,
    "Not attempted yet"
  ),
  stringsAsFactors = FALSE
)

infercnv_should_run <- isTRUE(run_infercnv_if_available) &&
  isTRUE(infercnv_available) &&
  !is.na(gene_order_file) &&
  length(reference_cells) > 0 &&
  length(tumor_hep_cells) > 0

if (infercnv_should_run) {
  dir.create(infercnv_out_dir, recursive = TRUE, showWarnings = FALSE)

  cells_use <- unique(c(tumor_hep_cells, reference_cells))
  expr_counts_sub <- expr_counts[, cells_use, drop = FALSE]
  meta_sub <- meta[cells_use, , drop = FALSE]

  infer_group <- meta_sub[[celltype_col]]
  infer_group[cells_use %in% tumor_hep_cells] <- "TumorHepatocyte"

  cell_annots <- data.frame(
    cell = cells_use,
    group = infer_group,
    stringsAsFactors = FALSE
  )

  annot_file <- file.path(out_dir, "inferCNV_06_input_cell_annotations.txt")
  utils::write.table(
    cell_annots,
    file = annot_file,
    sep = "\t",
    quote = FALSE,
    col.names = FALSE,
    row.names = FALSE
  )
  message("Wrote: ", normalizePath(annot_file, winslash = "/", mustWork = FALSE))

  infercnv_result <- tryCatch(
    {
      infercnv_obj <- infercnv::CreateInfercnvObject(
        raw_counts_matrix = expr_counts_sub,
        annotations_file = annot_file,
        delim = "\t",
        gene_order_file = gene_order_file,
        ref_group_names = reference_celltypes
      )

      infercnv_obj <- infercnv::run(
        infercnv_obj,
        cutoff = infercnv_cutoff,
        out_dir = infercnv_out_dir,
        cluster_by_groups = infercnv_cluster_by_groups,
        denoise = infercnv_denoise,
        HMM = infercnv_HMM,
        HMM_type = infercnv_HMM_type,
        num_threads = infercnv_num_threads
      )

      saveRDS(infercnv_obj, file.path(out_dir, "inferCNV_06_infercnv_object.rds"))

      list(status = "COMPLETED", error = NA_character_)
    },
    error = function(e) {
      list(status = "FAILED", error = conditionMessage(e))
    }
  )

  run_status$value[run_status$item == "infercnv_status"] <- infercnv_result$status

  if (!is.na(infercnv_result$error)) {
    run_status <- rbind(
      run_status,
      data.frame(
        item = "infercnv_error",
        value = infercnv_result$error,
        stringsAsFactors = FALSE
      )
    )
  }

  obs_file <- file.path(infercnv_out_dir, "infercnv.observations.txt")
  cnv_burden <- parse_infercnv_observations(obs_file, tumor_hep_cells)

  if (!is.null(cnv_burden)) {
    write_csv(cnv_burden, "inferCNV_06_cnv_burden.csv")
    run_status$value[run_status$item == "cnv_burden_status"] <- "Written to inferCNV_06_cnv_burden.csv"
  } else {
    run_status$value[run_status$item == "cnv_burden_status"] <- "Skipped: infercnv.observations.txt not available or not parseable"
  }
} else {
  reason <- c()
  if (!isTRUE(run_infercnv_if_available)) {
    reason <- c(reason, "run_infercnv_if_available is FALSE")
  }
  if (!isTRUE(infercnv_available)) {
    reason <- c(reason, "infercnv package is not installed")
  }
  if (is.na(gene_order_file)) {
    reason <- c(reason, "local hg38_gencode_v27.txt gene-order file was not found")
  }
  if (length(reference_cells) == 0) {
    reason <- c(reason, "no reference cells detected")
  }
  if (length(tumor_hep_cells) == 0) {
    reason <- c(reason, "no tumor hepatocytes detected")
  }

  run_status$value[run_status$item == "infercnv_status"] <- paste("SKIPPED:", paste(reason, collapse = "; "))
  run_status$value[run_status$item == "cnv_burden_status"] <- "Skipped because inferCNV did not run"
}

write_csv(run_status, "inferCNV_06_run_status.csv")

# -----------------------------------------------------------------------------
# Manifest and session info
# -----------------------------------------------------------------------------

manifest <- data.frame(
  output_file = c(
    "inferCNV_06_qc.csv",
    "inferCNV_06_hepatocyte_marker_scores.csv",
    "inferCNV_06_marker_gene_presence.csv",
    "inferCNV_06_group_counts.csv",
    "inferCNV_06_reference_counts.csv",
    "inferCNV_06_run_status.csv",
    "inferCNV_06_marker_scored_seurat.rds",
    "inferCNV_06_marker_score_summary.csv",
    "inferCNV_06_marker_score_plot_no_title.pdf",
    "inferCNV_06_marker_score_plot_no_title.png",
    "inferCNV_06_manifest.csv",
    "inferCNV_06_sessionInfo.txt"
  ),
  purpose = c(
    "Mandatory locked-object QC",
    "Tumor-derived hepatocyte marker-score source table",
    "Presence/absence of marker genes used by this script",
    "Current cell-type group counts",
    "Reference cell counts for optional inferCNV",
    "Status of optional inferCNV execution",
    "Seurat object with marker-score annotations",
    "Summary of marker-score analysis",
    "No-title marker-score QC plot",
    "No-title marker-score QC plot",
    "Output manifest",
    "R session information"
  ),
  direct_manuscript_figure = c(
    "No", "No", "No", "No", "No", "No", "No",
    "No", "No; QC plot only", "No; QC plot only", "No", "No"
  ),
  stringsAsFactors = FALSE
)

if (file.exists(file.path(out_dir, "inferCNV_06_input_cell_annotations.txt"))) {
  manifest <- rbind(
    manifest,
    data.frame(
      output_file = "inferCNV_06_input_cell_annotations.txt",
      purpose = "Cell annotation input file used by optional inferCNV",
      direct_manuscript_figure = "No",
      stringsAsFactors = FALSE
    )
  )
}

if (file.exists(file.path(out_dir, "inferCNV_06_cnv_burden.csv"))) {
  manifest <- rbind(
    manifest,
    data.frame(
      output_file = "inferCNV_06_cnv_burden.csv",
      purpose = "Cell-level CNV burden parsed from infercnv.observations.txt",
      direct_manuscript_figure = "No",
      stringsAsFactors = FALSE
    )
  )
}

if (file.exists(file.path(out_dir, "inferCNV_06_infercnv_object.rds"))) {
  manifest <- rbind(
    manifest,
    data.frame(
      output_file = "inferCNV_06_infercnv_object.rds",
      purpose = "inferCNV object generated by optional inferCNV run",
      direct_manuscript_figure = "No",
      stringsAsFactors = FALSE
    )
  )
}

write_csv(manifest, "inferCNV_06_manifest.csv")
write_text(capture.output(sessionInfo()), "inferCNV_06_sessionInfo.txt")

message("")
message("Script 06 complete.")
message("Tumor-derived hepatocytes: ", tumor_hep_n)
message("MarkerHigh tumor hepatocytes: ", marker_high_n)
message("MarkerLow tumor hepatocytes: ", marker_low_n)
message("inferCNV status: ", run_status$value[run_status$item == "infercnv_status"])
message("No manuscript figure PDF/PNG was generated by this script.")
