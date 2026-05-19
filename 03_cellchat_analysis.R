# =============================================================================
# 03_cellchat_analysis.R
#
# Current purpose:
#   Current-object CellChat analysis with locked metadata/QC guard.
#
# This replacement removes legacy local paths and obsolete labels:
#   - no setwd("D:/scRNA_project")
#   - no tumor_hepatocytes.rds dependency
#   - no "Hepatocytes" plural label for tumor-derived hepatocytes
#
# Locked current-object convention:
#   input object      = seurat_final.rds
#   patient column    = patient
#   cell-type column  = cell_type
#   site column       = site
#   glycolysis column = Glycolysis_AUC
#
# Tumor-derived hepatocytes:
#   site == "Tumor" & cell_type == "Hepatocyte"
#   n = 15,391
#   median Glycolysis_AUC cutoff = 0.2203849
#   GlycoHigh = Glycolysis_AUC > median_cut, n = 7,695
#   GlycoLow  = remaining tumor-derived hepatocytes, n = 7,696
#
# CellChat rounds:
#   Round 1: original cell_type communication landscape
#   Round 2: GlycoHigh/GlycoLow tumor hepatocytes + remaining cell types
#   Round 3: optional subtype-resolved run only if an explicit subtype column
#            already exists in the object. This script does not de novo
#            recluster T/NK or myeloid cells because that would create another
#            moving target unless a subtype annotation has been QC-locked.
#
# Outputs:
#   results/CellChat_03_qc.csv
#   results/CellChat_03_input_cell_metadata.csv
#   results/CellChat_03_round1_group_counts.csv
#   results/CellChat_03_round2_group_counts.csv
#   results/CellChat_03_round1.rds
#   results/CellChat_03_round2_glyco.rds
#   results/CellChat_03_round1_communication.csv
#   results/CellChat_03_round2_glyco_communication.csv
#   results/CellChat_03_round1_net_count.csv
#   results/CellChat_03_round1_net_weight.csv
#   results/CellChat_03_round2_glyco_net_count.csv
#   results/CellChat_03_round2_glyco_net_weight.csv
#   results/CellChat_03_round3_status.csv
#   results/CellChat_03_manifest.csv
#   results/CellChat_03_sessionInfo.txt
#
# Notes:
#   This script does not generate manuscript figures directly.
#   Figure generation should use exported source/QC tables or locked
#   downstream figure-specific scripts.
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

input_rds <- "seurat_final.rds"
out_dir <- "results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

patient_col <- "patient"
celltype_col <- "cell_type"
site_col <- "site"
gly_col <- "Glycolysis_AUC"

expected_tumor_hepatocytes <- 15391L
expected_median_cut <- 0.2203849
expected_glyco_high <- 7695L
expected_glyco_low <- 7696L

cellchat_min_cells <- 10L
cellchat_nboot <- 1000L

# If TRUE, restrict CellChat database to secreted signaling interactions.
# Keep FALSE to preserve the broader original CellChat landscape behavior.
use_secreted_signaling_only <- FALSE

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

matrix_to_long <- function(mat, value_col = "value") {
  if (is.null(mat)) {
    return(data.frame(source = character(), target = character(), value = numeric()))
  }
  df <- as.data.frame(as.table(as.matrix(mat)), stringsAsFactors = FALSE)
  names(df) <- c("source", "target", value_col)
  df
}

required_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Required package is not installed: ", pkg, call. = FALSE)
  }
}

get_expression_matrix <- function(seu, assay_use = NULL) {
  if (is.null(assay_use)) {
    assay_use <- Seurat::DefaultAssay(seu)
    if ("RNA" %in% Seurat::Assays(seu)) {
      assay_use <- "RNA"
    }
  }

  message("Using assay: ", assay_use)

  expr <- NULL

  # SeuratObject::LayerData is preferred for Seurat v5 when available.
  expr <- tryCatch(
    {
      if (requireNamespace("SeuratObject", quietly = TRUE) &&
          "LayerData" %in% getNamespaceExports("SeuratObject")) {
        SeuratObject::LayerData(seu, assay = assay_use, layer = "data")
      } else {
        NULL
      }
    },
    error = function(e) NULL
  )

  # Seurat v5 GetAssayData(layer = "data")
  if (is.null(expr)) {
    expr <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay_use, layer = "data"),
      error = function(e) NULL
    )
  }

  # Seurat v4 fallback
  if (is.null(expr)) {
    expr <- tryCatch(
      Seurat::GetAssayData(seu, assay = assay_use, slot = "data"),
      error = function(e) NULL
    )
  }

  if (is.null(expr)) {
    stop("Could not retrieve normalized expression matrix from assay: ", assay_use)
  }

  expr
}

make_group_counts <- function(meta, group_col) {
  tab <- sort(table(meta[[group_col]]), decreasing = TRUE)
  data.frame(
    group = names(tab),
    n_cells = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

run_cellchat_round <- function(expr, meta, group_col, prefix,
                               min_cells = 10L,
                               nboot = 1000L,
                               use_secreted_only = FALSE) {
  message("")
  message("Running CellChat: ", prefix)
  message("Grouping column: ", group_col)

  group_counts <- make_group_counts(meta, group_col)
  write_csv(group_counts, paste0(prefix, "_group_counts.csv"))

  if (nrow(group_counts) < 2) {
    stop("CellChat requires at least two groups for ", prefix, call. = FALSE)
  }

  cellchat <- CellChat::createCellChat(
    object = expr,
    meta = meta,
    group.by = group_col
  )

  db_use <- CellChat::CellChatDB.human
  if (isTRUE(use_secreted_only)) {
    db_use <- CellChat::subsetDB(db_use, search = "Secreted Signaling")
  }
  cellchat@DB <- db_use

  cellchat <- CellChat::subsetData(cellchat)
  cellchat <- CellChat::identifyOverExpressedGenes(cellchat)
  cellchat <- CellChat::identifyOverExpressedInteractions(cellchat)

  cellchat <- CellChat::computeCommunProb(
    cellchat,
    type = "triMean",
    nboot = nboot
  )

  cellchat <- CellChat::filterCommunication(
    cellchat,
    min.cells = min_cells
  )

  cellchat <- CellChat::computeCommunProbPathway(cellchat)
  cellchat <- CellChat::aggregateNet(cellchat)

  rds_file <- file.path(out_dir, paste0(prefix, ".rds"))
  saveRDS(cellchat, rds_file)
  message("Wrote: ", normalizePath(rds_file, winslash = "/", mustWork = FALSE))

  comm_df <- tryCatch(
    CellChat::subsetCommunication(cellchat),
    error = function(e) {
      warning("subsetCommunication failed for ", prefix, ": ", conditionMessage(e))
      data.frame()
    }
  )
  write_csv(comm_df, paste0(prefix, "_communication.csv"))

  count_long <- matrix_to_long(cellchat@net$count, value_col = "n_interactions")
  weight_long <- matrix_to_long(cellchat@net$weight, value_col = "interaction_weight")

  write_csv(count_long, paste0(prefix, "_net_count.csv"))
  write_csv(weight_long, paste0(prefix, "_net_weight.csv"))

  invisible(list(
    object = cellchat,
    rds_file = rds_file,
    communication_n = nrow(comm_df),
    group_counts = group_counts
  ))
}

# -----------------------------------------------------------------------------
# Package checks
# -----------------------------------------------------------------------------

required_package("Seurat")
required_package("CellChat")

suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
})

# -----------------------------------------------------------------------------
# Load and validate input object
# -----------------------------------------------------------------------------

if (!file.exists(input_rds)) {
  stop("Input object not found: ", input_rds, call. = FALSE)
}

seu <- readRDS(input_rds)
meta <- seu@meta.data

required_cols <- c(patient_col, celltype_col, site_col, gly_col)
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

if (!is.numeric(meta[[gly_col]])) {
  meta[[gly_col]] <- suppressWarnings(as.numeric(meta[[gly_col]]))
}

if (anyNA(meta[[gly_col]])) {
  warning("Some Glycolysis_AUC values are NA after numeric conversion.")
}

tumor_hepa <- meta[[site_col]] == "Tumor" &
  meta[[celltype_col]] == "Hepatocyte" &
  !is.na(meta[[gly_col]])

tumor_hepa_n <- sum(tumor_hepa)
median_cut <- stats::median(meta[[gly_col]][tumor_hepa], na.rm = TRUE)

glyco_group <- rep(NA_character_, nrow(meta))
glyco_group[tumor_hepa & meta[[gly_col]] > median_cut] <- "TumorGlycoHigh"
glyco_group[tumor_hepa & meta[[gly_col]] <= median_cut] <- "TumorGlycoLow"

cell_type_glyco <- meta[[celltype_col]]
cell_type_glyco[tumor_hepa] <- glyco_group[tumor_hepa]

meta$cell_type_current <- meta[[celltype_col]]
meta$cell_type_glyco <- safe_chr(cell_type_glyco)
meta$tumor_hepatocyte_status <- ifelse(tumor_hepa, "TumorHepatocyte", "OtherCell")
meta$glyco_group_tumor_hepatocyte <- ifelse(tumor_hepa, glyco_group, NA_character_)

glyco_high_n <- sum(meta$cell_type_glyco == "TumorGlycoHigh")
glyco_low_n <- sum(meta$cell_type_glyco == "TumorGlycoLow")

qc_df <- data.frame(
  check_id = c(
    "input_rds",
    "patient_column_present",
    "celltype_column_present",
    "site_column_present",
    "glycolysis_column_present",
    "tumor_hepatocyte_definition",
    "tumor_hepatocyte_n",
    "median_cutoff",
    "glycohigh_definition",
    "glycolow_definition",
    "glycohigh_n",
    "glycolow_n",
    "legacy_plural_hepatocytes_not_used",
    "no_external_tumor_hepatocytes_rds_dependency"
  ),
  expected = c(
    "seurat_final.rds",
    patient_col,
    celltype_col,
    site_col,
    gly_col,
    'site == "Tumor" & cell_type == "Hepatocyte"',
    as.character(expected_tumor_hepatocytes),
    sprintf("%.7f", expected_median_cut),
    "Glycolysis_AUC > median_cut",
    "remaining tumor-derived hepatocytes",
    as.character(expected_glyco_high),
    as.character(expected_glyco_low),
    'cell_type == "Hepatocyte"',
    "No"
  ),
  observed = c(
    input_rds,
    ifelse(patient_col %in% colnames(meta), patient_col, "MISSING"),
    ifelse(celltype_col %in% colnames(meta), celltype_col, "MISSING"),
    ifelse(site_col %in% colnames(meta), site_col, "MISSING"),
    ifelse(gly_col %in% colnames(meta), gly_col, "MISSING"),
    'site == "Tumor" & cell_type == "Hepatocyte"',
    as.character(tumor_hepa_n),
    sprintf("%.7f", median_cut),
    "Glycolysis_AUC > median_cut",
    "remaining tumor-derived hepatocytes",
    as.character(glyco_high_n),
    as.character(glyco_low_n),
    'cell_type == "Hepatocyte"',
    "No"
  ),
  pass = c(
    identical(input_rds, "seurat_final.rds"),
    patient_col %in% colnames(meta),
    celltype_col %in% colnames(meta),
    site_col %in% colnames(meta),
    gly_col %in% colnames(meta),
    TRUE,
    identical(as.integer(tumor_hepa_n), expected_tumor_hepatocytes),
    isTRUE(abs(median_cut - expected_median_cut) < 1e-6),
    TRUE,
    TRUE,
    identical(as.integer(glyco_high_n), expected_glyco_high),
    identical(as.integer(glyco_low_n), expected_glyco_low),
    !any(meta[[celltype_col]][tumor_hepa] == "Hepatocytes"),
    TRUE
  ),
  stringsAsFactors = FALSE
)

write_csv(qc_df, "CellChat_03_qc.csv")

if (!all(qc_df$pass)) {
  stop(
    "Mandatory CellChat input QC failed. See results/CellChat_03_qc.csv. ",
    "Do not run CellChat on a non-locked object/annotation.",
    call. = FALSE
  )
}

# -----------------------------------------------------------------------------
# Build CellChat input metadata and expression matrix
# -----------------------------------------------------------------------------

expr <- get_expression_matrix(seu)

common_cells <- intersect(colnames(expr), rownames(meta))
if (length(common_cells) == 0) {
  stop("No overlapping cell barcodes between expression matrix and metadata.")
}

expr <- expr[, common_cells, drop = FALSE]
meta <- meta[common_cells, , drop = FALSE]

input_meta <- data.frame(
  cell = rownames(meta),
  patient = meta[[patient_col]],
  site = meta[[site_col]],
  cell_type = meta[[celltype_col]],
  Glycolysis_AUC = meta[[gly_col]],
  tumor_hepatocyte_status = meta$tumor_hepatocyte_status,
  glyco_group_tumor_hepatocyte = meta$glyco_group_tumor_hepatocyte,
  cell_type_current = meta$cell_type_current,
  cell_type_glyco = meta$cell_type_glyco,
  stringsAsFactors = FALSE
)

write_csv(input_meta, "CellChat_03_input_cell_metadata.csv")

write_csv(make_group_counts(meta, "cell_type_current"), "CellChat_03_round1_group_counts.csv")
write_csv(make_group_counts(meta, "cell_type_glyco"), "CellChat_03_round2_group_counts.csv")

# -----------------------------------------------------------------------------
# Run CellChat: Round 1 and Round 2
# -----------------------------------------------------------------------------

round1 <- run_cellchat_round(
  expr = expr,
  meta = meta,
  group_col = "cell_type_current",
  prefix = "CellChat_03_round1",
  min_cells = cellchat_min_cells,
  nboot = cellchat_nboot,
  use_secreted_only = use_secreted_signaling_only
)

round2 <- run_cellchat_round(
  expr = expr,
  meta = meta,
  group_col = "cell_type_glyco",
  prefix = "CellChat_03_round2_glyco",
  min_cells = cellchat_min_cells,
  nboot = cellchat_nboot,
  use_secreted_only = use_secreted_signaling_only
)

# -----------------------------------------------------------------------------
# Optional Round 3: subtype-resolved only if an explicit subtype column exists
# -----------------------------------------------------------------------------

subtype_candidates <- c(
  "cell_subtype",
  "subtype",
  "Cell_subtype",
  "immune_subtype",
  "major_subtype",
  "annotation_subtype"
)

subtype_col <- subtype_candidates[subtype_candidates %in% colnames(meta)][1]

round3_status <- data.frame(
  item = c(
    "round3_requested",
    "subtype_column_detected",
    "subtype_column_used",
    "de_novo_subclustering",
    "status"
  ),
  value = c(
    "optional",
    ifelse(is.na(subtype_col), "No", "Yes"),
    ifelse(is.na(subtype_col), "None", subtype_col),
    "No",
    ifelse(
      is.na(subtype_col),
      "Skipped: no explicit QC-locked subtype column found in seurat_final.rds.",
      "Will run using the detected explicit subtype column."
    )
  ),
  stringsAsFactors = FALSE
)

write_csv(round3_status, "CellChat_03_round3_status.csv")

round3 <- NULL
if (!is.na(subtype_col)) {
  meta$cell_subtype_glyco <- meta$cell_type_glyco

  subtype_values <- safe_chr(meta[[subtype_col]])
  use_subtype <- meta[[celltype_col]] %in% c(
    "T_NK", "T/NK", "CD8T", "CD8_T", "T", "NK",
    "Myeloid", "Macrophage", "Monocyte", "DC"
  )

  meta$cell_subtype_glyco[use_subtype] <- subtype_values[use_subtype]
  meta$cell_subtype_glyco <- safe_chr(meta$cell_subtype_glyco)

  write_csv(make_group_counts(meta, "cell_subtype_glyco"),
            "CellChat_03_round3_group_counts.csv")

  round3 <- run_cellchat_round(
    expr = expr,
    meta = meta,
    group_col = "cell_subtype_glyco",
    prefix = "CellChat_03_round3_subtype_glyco",
    min_cells = cellchat_min_cells,
    nboot = cellchat_nboot,
    use_secreted_only = use_secreted_signaling_only
  )
}

# -----------------------------------------------------------------------------
# Manifest and session info
# -----------------------------------------------------------------------------

manifest <- data.frame(
  output_file = c(
    "CellChat_03_qc.csv",
    "CellChat_03_input_cell_metadata.csv",
    "CellChat_03_round1_group_counts.csv",
    "CellChat_03_round2_group_counts.csv",
    "CellChat_03_round1.rds",
    "CellChat_03_round2_glyco.rds",
    "CellChat_03_round1_communication.csv",
    "CellChat_03_round2_glyco_communication.csv",
    "CellChat_03_round1_net_count.csv",
    "CellChat_03_round1_net_weight.csv",
    "CellChat_03_round2_glyco_net_count.csv",
    "CellChat_03_round2_glyco_net_weight.csv",
    "CellChat_03_round3_status.csv",
    "CellChat_03_sessionInfo.txt"
  ),
  purpose = c(
    "Mandatory locked-object QC",
    "Cell-level metadata source table used for CellChat",
    "Round 1 group counts",
    "Round 2 GlycoHigh/GlycoLow group counts",
    "Round 1 CellChat object",
    "Round 2 CellChat object",
    "Round 1 inferred ligand-receptor communication table",
    "Round 2 inferred ligand-receptor communication table",
    "Round 1 source-target interaction-count matrix in long format",
    "Round 1 source-target interaction-weight matrix in long format",
    "Round 2 source-target interaction-count matrix in long format",
    "Round 2 source-target interaction-weight matrix in long format",
    "Round 3 subtype-resolution status",
    "R session information"
  ),
  figure_file = "No direct manuscript figure is generated by this script.",
  stringsAsFactors = FALSE
)

if (!is.null(round3)) {
  manifest <- rbind(
    manifest,
    data.frame(
      output_file = c(
        "CellChat_03_round3_group_counts.csv",
        "CellChat_03_round3_subtype_glyco.rds",
        "CellChat_03_round3_subtype_glyco_communication.csv",
        "CellChat_03_round3_subtype_glyco_net_count.csv",
        "CellChat_03_round3_subtype_glyco_net_weight.csv"
      ),
      purpose = c(
        "Round 3 subtype-resolved group counts",
        "Round 3 subtype-resolved CellChat object",
        "Round 3 inferred ligand-receptor communication table",
        "Round 3 source-target interaction-count matrix in long format",
        "Round 3 source-target interaction-weight matrix in long format"
      ),
      figure_file = "No direct manuscript figure is generated by this script.",
      stringsAsFactors = FALSE
    )
  )
}

write_csv(manifest, "CellChat_03_manifest.csv")

session_lines <- capture.output(sessionInfo())
write_text(session_lines, "CellChat_03_sessionInfo.txt")

message("")
message("03 CellChat analysis complete.")
message("Mandatory locked-object QC passed.")
message("Round 1 communication rows: ", round1$communication_n)
message("Round 2 communication rows: ", round2$communication_n)
if (is.null(round3)) {
  message("Round 3 skipped because no explicit QC-locked subtype column was detected.")
} else {
  message("Round 3 communication rows: ", round3$communication_n)
}
message("No manuscript figure PDF/PNG was generated by this script.")
