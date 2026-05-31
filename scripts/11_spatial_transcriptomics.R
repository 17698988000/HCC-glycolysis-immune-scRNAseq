# ============================================================
# 11_spatial_transcriptomics.R
#
# Purpose:
#   Spatial transcriptomics validation of glycolysis-linked
#   immunosuppressive ligand programs in HCC Visium samples.
#
# Final manuscript alignment:
#   Dataset: GSE238264
#   Samples: HCC1R, HCC2R, HCC3R, HCC4R
#   Method:
#     - Load 10x Visium data
#     - QC: nFeature_Spatial >= 200 and percent.mt < 25
#     - Normalize with SCTransform
#     - Compute Glycolysis1 module score using the 22-gene glycolysis set
#     - Split each sample into spatial GlycoHigh/GlycoLow by within-sample median Glycolysis1
#     - Evaluate ENO1, four-gene glycolysis score, PTGES, MIF, SPP1,
#       and MIF/SPP1 ligand score
#
# Important interpretation:
#   Visium spots are mixed-cell tissue neighborhoods.
#   These analyses support spot-level tissue co-enrichment only.
#   They do not prove same-cell co-expression, tumor-cell-specific ligand
#   production, direct ligand-receptor contact, or causal signaling.
#
# Output policy:
#   1. Write source CSV.
#   2. Write summary CSV.
#   3. Write QC CSV.
#   4. Generate PDF/PNG only if mandatory QC checks pass.
#
# ============================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0. Package checks
# -----------------------------

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "readr",
  "ggplot2",
  "patchwork",
  "purrr"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall missing CRAN/Bioconductor packages before running this script."
  )
}

library(Seurat)
library(SeuratObject)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(patchwork)
library(purrr)

# Avoid common namespace masking
select      <- dplyr::select
mutate      <- dplyr::mutate
filter      <- dplyr::filter
arrange     <- dplyr::arrange
summarise   <- dplyr::summarise
group_by    <- dplyr::group_by
ungroup     <- dplyr::ungroup
left_join   <- dplyr::left_join
bind_rows   <- dplyr::bind_rows
all_of      <- dplyr::all_of

set.seed(20260519)

# -----------------------------
# 1. User parameters
# -----------------------------

data_dir <- file.path("data", "GSE238264")
outdir <- file.path("results", "GSE238264_spatial")

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sample_ids <- c("HCC1R", "HCC2R", "HCC3R", "HCC4R")

# Representative sample for Figure 9A spatial feature maps.
representative_sample <- "HCC1R"

# QC cutoffs from final Methods.
min_features <- 200
max_percent_mt <- 25

glyco_genes <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKM", "ALDOA", "ALDOB", "ALDOC",
  "TPI1", "GAPDH", "PGK1", "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
  "SLC2A1", "SLC2A3", "PFKFB3", "GCK"
)

four_gene_coef <- c(
  TPI1   = 0.3041908,
  ENO1   = 0.9639654,
  LDHA   = 1.3404374,
  SLC2A1 = 0.2424239
)

feature_genes <- c("ENO1", "LDHA", "SPP1", "MIF", "PTGES", names(four_gene_coef))

outcome_features <- c(
  "ENO1",
  "four_gene_score",
  "PTGES",
  "MIF",
  "MIF_SPP1_ligand_score",
  "SPP1"
)

expected_group_counts <- tibble::tribble(
  ~sample, ~expected_high, ~expected_low,
  "HCC1R", 1503L, 1503L,
  "HCC2R", 1383L, 1383L,
  "HCC3R", 1085L, 1085L,
  "HCC4R", 1501L, 1501L
)

expected_median_effects <- tibble::tribble(
  ~feature, ~expected_median_effect,
  "ENO1", 0.911,
  "four_gene_score", 0.911,
  "PTGES", 0.444,
  "MIF", 0.347,
  "MIF_SPP1_ligand_score", 0.329,
  "SPP1", 0.261
)

# Optional RCTD composition-adjusted regression hook.
# If this file exists, it should contain one row per spot with columns:
# sample, cell or barcode, and broad cell-type proportions:
# Hepatocyte/hepatocyte, Myeloid/myeloid, T_NK/TNK/T.NK, Fibroblast/fibroblast,
# Endothelial/endothelial, B_cell/B.cells/B.
rctd_proportion_file <- file.path(data_dir, "rctd_spot_proportions.csv")

# -----------------------------
# 2. Helper functions
# -----------------------------

qc_row <- function(check, observed, expected, pass, severity = "ERROR", notes = "") {
  tibble(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = as.logical(pass),
    severity = severity,
    notes = notes
  )
}

get_assay_data_safe <- function(seu, assay, layer = "data", slot = "data") {
  tryCatch(
    {
      SeuratObject::LayerData(seu, assay = assay, layer = layer)
    },
    error = function(e) {
      Seurat::GetAssayData(seu, assay = assay, slot = slot)
    }
  )
}

file_exists_any <- function(paths) {
  any(file.exists(paths))
}

has_matrix_tsv <- function(path) {
  file_exists_any(file.path(path, c("matrix.mtx", "matrix.mtx.gz"))) &&
    file_exists_any(file.path(path, c("features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"))) &&
    file_exists_any(file.path(path, c("barcodes.tsv", "barcodes.tsv.gz")))
}

is_visium_dir <- function(path) {
  if (!dir.exists(path)) return(FALSE)

  has_spatial <- dir.exists(file.path(path, "spatial"))
  has_h5 <- file.exists(file.path(path, "filtered_feature_bc_matrix.h5"))
  has_mtx_root <- has_matrix_tsv(path)
  has_mtx_filtered <- has_matrix_tsv(file.path(path, "filtered_feature_bc_matrix"))
  has_mtx_raw <- has_matrix_tsv(file.path(path, "raw_feature_bc_matrix"))

  has_spatial && (has_h5 || has_mtx_root || has_mtx_filtered || has_mtx_raw)
}

find_sample_dir <- function(data_dir, sid) {
  candidates <- unique(c(
    file.path(data_dir, sid),
    file.path(data_dir, sid, sid),
    file.path(data_dir, paste0(sid, "_outs")),
    file.path(data_dir, paste0(sid, "_filtered")),
    data_dir
  ))

  valid <- candidates[vapply(candidates, is_visium_dir, logical(1))]

  if (length(valid) == 0) {
    stop(
      "Could not locate a valid Visium directory for sample ", sid, ".\n",
      "Expected one of these layouts:\n",
      "  data/GSE238264/", sid, "/filtered_feature_bc_matrix.h5 plus spatial/\n",
      "  data/GSE238264/", sid, "/filtered_feature_bc_matrix/ plus spatial/\n",
      "  data/GSE238264/", sid, "/", sid, "/filtered_feature_bc_matrix.h5 plus spatial/\n",
      "Checked candidates:\n",
      paste(candidates, collapse = "\n")
    )
  }

  valid[1]
}

load_visium_sample <- function(sample_dir, sid) {
  cat("\nLoading sample ", sid, " from: ", sample_dir, "\n", sep = "")

  # Primary Seurat loader for standard Space Ranger output.
  seu <- tryCatch(
    {
      Seurat::Load10X_Spatial(
        data.dir = sample_dir,
        slice = sid,
        assay = "Spatial",
        filter.matrix = TRUE
      )
    },
    error = function(e) {
      message("Load10X_Spatial failed for ", sid, ": ", conditionMessage(e))
      NULL
    }
  )

  if (!is.null(seu)) {
    return(seu)
  }

  # Fallback for matrix.mtx-style folders.
  matrix_dir <- NULL

  if (has_matrix_tsv(file.path(sample_dir, "filtered_feature_bc_matrix"))) {
    matrix_dir <- file.path(sample_dir, "filtered_feature_bc_matrix")
  } else if (has_matrix_tsv(file.path(sample_dir, "raw_feature_bc_matrix"))) {
    matrix_dir <- file.path(sample_dir, "raw_feature_bc_matrix")
  } else if (has_matrix_tsv(sample_dir)) {
    matrix_dir <- sample_dir
  }

  if (is.null(matrix_dir)) {
    stop("Could not find readable matrix.mtx/features.tsv/barcodes.tsv files for ", sid)
  }

  counts <- Seurat::Read10X(data.dir = matrix_dir)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  seu <- Seurat::CreateSeuratObject(
    counts = counts,
    assay = "Spatial",
    project = sid
  )

  warning(
    "Loaded ", sid, " with Read10X fallback. Spatial images may not be available; ",
    "source/QC tables will still be generated, but SpatialFeaturePlot may be replaced by coordinate plots if coordinates are readable."
  )

  spatial_dir <- file.path(sample_dir, "spatial")
  coord_candidates <- file.path(
    spatial_dir,
    c("tissue_positions.csv", "tissue_positions_list.csv")
  )
  coord_file <- coord_candidates[file.exists(coord_candidates)][1]

  if (!is.na(coord_file)) {
    coords <- readr::read_csv(coord_file, col_names = FALSE, show_col_types = FALSE)

    if (ncol(coords) >= 6) {
      names(coords)[1:6] <- c("cell", "in_tissue", "array_row", "array_col", "imagerow", "imagecol")
      coords <- coords %>%
        filter(cell %in% colnames(seu)) %>%
        select(cell, in_tissue, array_row, array_col, imagerow, imagecol)

      rownames(coords) <- coords$cell
      seu <- AddMetaData(seu, metadata = coords[colnames(seu), c("in_tissue", "array_row", "array_col", "imagerow", "imagecol")])
    }
  }

  seu
}

get_spatial_coordinates_safe <- function(seu) {
  coords <- tryCatch(
    {
      Seurat::GetTissueCoordinates(seu) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("cell")
    },
    error = function(e) {
      NULL
    }
  )

  if (!is.null(coords) && nrow(coords) > 0) {
    return(coords)
  }

  meta <- seu@meta.data %>%
    tibble::rownames_to_column("cell")

  coord_cols <- intersect(c("imagerow", "imagecol", "array_row", "array_col"), names(meta))

  if (length(coord_cols) >= 2) {
    return(meta %>% select(cell, all_of(coord_cols)))
  }

  tibble(cell = colnames(seu))
}

extract_expr_df <- function(seu, sample_id, genes, assay = "SCT") {
  mat <- get_assay_data_safe(seu, assay = assay, layer = "data", slot = "data")

  present <- intersect(unique(genes), rownames(mat))
  missing <- setdiff(unique(genes), rownames(mat))

  if (length(missing) > 0) {
    warning("Missing genes in ", sample_id, ": ", paste(missing, collapse = ", "))
  }

  expr <- as.matrix(mat[present, , drop = FALSE]) %>%
    t() %>%
    as.data.frame(check.names = FALSE) %>%
    tibble::rownames_to_column("cell")

  for (g in unique(genes)) {
    if (!g %in% names(expr)) {
      expr[[g]] <- NA_real_
    }
  }

  expr %>% select(cell, all_of(unique(genes)))
}

make_coordinate_feature_plot <- function(df, feature, title_text) {
  if (all(c("imagecol", "imagerow") %in% names(df))) {
    x_col <- "imagecol"
    y_col <- "imagerow"
  } else if (all(c("x", "y") %in% names(df))) {
    x_col <- "x"
    y_col <- "y"
  } else if (all(c("array_col", "array_row") %in% names(df))) {
    x_col <- "array_col"
    y_col <- "array_row"
  } else {
    stop("No coordinate columns available for fallback coordinate plotting.")
  }

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = .data[[feature]])) +
    geom_point(size = 0.45) +
    scale_y_reverse() +
    coord_fixed() +
    labs(title = title_text, x = NULL, y = NULL, color = feature) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

make_coordinate_group_plot <- function(df, title_text) {
  if (all(c("imagecol", "imagerow") %in% names(df))) {
    x_col <- "imagecol"
    y_col <- "imagerow"
  } else if (all(c("x", "y") %in% names(df))) {
    x_col <- "x"
    y_col <- "y"
  } else if (all(c("array_col", "array_row") %in% names(df))) {
    x_col <- "array_col"
    y_col <- "array_row"
  } else {
    stop("No coordinate columns available for fallback coordinate plotting.")
  }

  ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], color = GlycoGroup)) +
    geom_point(size = 0.45, alpha = 0.9) +
    scale_y_reverse() +
    coord_fixed() +
    scale_color_manual(values = c("GlycoLow" = "#4575B4", "GlycoHigh" = "#D73027")) +
    labs(title = title_text, x = NULL, y = NULL, color = NULL) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
      legend.position = "right"
    )
}

standardize_rctd_cols <- function(df) {
  cn <- names(df)
  cn_clean <- tolower(stringr::str_replace_all(cn, "[^a-z0-9]+", "_"))

  names(df) <- cn_clean

  sample_col <- intersect(c("sample", "sample_id", "sid"), names(df))[1]
  cell_col <- intersect(c("cell", "barcode", "spot", "spot_barcode"), names(df))[1]

  if (is.na(sample_col) || is.na(cell_col)) {
    stop("RCTD proportion file must contain sample and cell/barcode columns.")
  }

  df <- df %>%
    rename(sample = all_of(sample_col), cell = all_of(cell_col))

  map_one <- function(candidates, target) {
    hit <- intersect(candidates, names(df))[1]
    if (!is.na(hit)) {
      df[[target]] <<- df[[hit]]
    } else {
      df[[target]] <<- NA_real_
    }
  }

  map_one(c("hepatocyte", "hepatocytes"), "prop_hepatocyte")
  map_one(c("myeloid", "myeloid_cells", "macrophage", "macrophages"), "prop_myeloid")
  map_one(c("t_nk", "tnk", "t_nk_cells", "t_cells_nk_cells"), "prop_t_nk")
  map_one(c("fibroblast", "fibroblasts"), "prop_fibroblast")
  map_one(c("endothelial", "endothelial_cells"), "prop_endothelial")
  map_one(c("b_cell", "b_cells", "b"), "prop_b_cell")

  df %>%
    select(
      sample, cell,
      prop_hepatocyte, prop_myeloid, prop_t_nk,
      prop_fibroblast, prop_endothelial, prop_b_cell
    )
}

run_rctd_regression_if_available <- function(source_df, rctd_file, outdir) {
  if (!file.exists(rctd_file)) {
    message("Optional RCTD proportion file not found: ", rctd_file)
    message("Skipping composition-adjusted regression. This does not affect primary spatial QC.")
    return(NULL)
  }

  rctd <- readr::read_csv(rctd_file, show_col_types = FALSE) %>%
    standardize_rctd_cols()

  merged <- source_df %>%
    left_join(rctd, by = c("sample", "cell")) %>%
    filter(
      !is.na(prop_hepatocyte),
      !is.na(prop_myeloid),
      !is.na(prop_t_nk),
      !is.na(prop_fibroblast),
      !is.na(prop_endothelial)
    ) %>%
    group_by(sample) %>%
    mutate(
      Glycolysis1_z = as.numeric(scale(Glycolysis1)),
      log_depth = log1p(nCount_Spatial)
    ) %>%
    ungroup()

  if (nrow(merged) == 0) {
    warning("RCTD file was found, but no rows could be merged with spatial source data.")
    return(NULL)
  }

  results <- bind_rows(lapply(outcome_features, function(feat) {
    dat <- merged %>%
      filter(!is.na(.data[[feat]])) %>%
      group_by(sample) %>%
      mutate(outcome_z = as.numeric(scale(.data[[feat]]))) %>%
      ungroup() %>%
      filter(!is.na(outcome_z), !is.na(Glycolysis1_z))

    if (nrow(dat) < 20) {
      return(tibble(
        feature = feat,
        model = "combined_RCTD_adjusted",
        n_spots = nrow(dat),
        beta = NA_real_,
        p_value = NA_real_,
        notes = "Insufficient merged spots"
      ))
    }

    fit <- lm(
      outcome_z ~ Glycolysis1_z +
        prop_hepatocyte +
        prop_myeloid +
        prop_t_nk +
        prop_fibroblast +
        prop_endothelial +
        log_depth +
        factor(sample),
      data = dat
    )

    coef_tab <- summary(fit)$coefficients

    tibble(
      feature = feat,
      model = "combined_RCTD_adjusted",
      n_spots = nrow(dat),
      beta = unname(coef_tab["Glycolysis1_z", "Estimate"]),
      p_value = unname(coef_tab["Glycolysis1_z", "Pr(>|t|)"]),
      notes = ""
    )
  })) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"))

  out_file <- file.path(outdir, "Fig9_spatial_RCTD_adjusted_regression_optional.csv")
  readr::write_csv(results, out_file)

  message("Optional RCTD-adjusted regression written to: ", out_file)

  results
}

# -----------------------------
# 3. Process Visium samples
# -----------------------------

spot_source_list <- list()
summary_list <- list()
feature_plot_objects <- list()
group_plot_objects <- list()

for (sid in sample_ids) {
  sample_dir <- find_sample_dir(data_dir, sid)
  seu <- load_visium_sample(sample_dir, sid)

  if (!"Spatial" %in% names(seu@assays)) {
    stop("Sample ", sid, " does not contain a Spatial assay.")
  }

  DefaultAssay(seu) <- "Spatial"

  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    assay = "Spatial",
    pattern = "^MT-"
  )

  pre_qc_n <- ncol(seu)

  seu <- subset(
    seu,
    subset = nFeature_Spatial >= min_features & percent.mt < max_percent_mt
  )

  post_qc_n <- ncol(seu)

  cat("Sample ", sid, ": pre-QC spots = ", pre_qc_n, "; post-QC spots = ", post_qc_n, "\n", sep = "")

  seu <- SCTransform(
    seu,
    assay = "Spatial",
    verbose = FALSE
  )

  DefaultAssay(seu) <- "SCT"

  glyco_present <- intersect(glyco_genes, rownames(seu))
  glyco_missing <- setdiff(glyco_genes, rownames(seu))

  if (length(glyco_present) < 10) {
    stop(
      "Too few glycolysis genes detected in ", sid, ": ",
      length(glyco_present), " present."
    )
  }

  seu <- AddModuleScore(
    seu,
    features = list(glyco_present),
    name = "Glycolysis",
    assay = "SCT",
    ctrl = 100
  )

  if (!"Glycolysis1" %in% colnames(seu@meta.data)) {
    stop("AddModuleScore did not create Glycolysis1 for ", sid)
  }

  med_score <- median(seu$Glycolysis1, na.rm = TRUE)

  seu$GlycoGroup <- factor(
    ifelse(seu$Glycolysis1 >= med_score, "GlycoHigh", "GlycoLow"),
    levels = c("GlycoLow", "GlycoHigh")
  )

  expr_df <- extract_expr_df(
    seu = seu,
    sample_id = sid,
    genes = unique(feature_genes),
    assay = "SCT"
  )

  coords <- get_spatial_coordinates_safe(seu)

  meta_df <- seu@meta.data %>%
    tibble::rownames_to_column("cell") %>%
    mutate(
      sample = sid,
      pre_qc_spots = pre_qc_n,
      post_qc_spots = post_qc_n,
      median_Glycolysis1 = med_score,
      n_glyco_genes_present = length(glyco_present),
      glyco_genes_missing = paste(glyco_missing, collapse = ";")
    ) %>%
    select(
      sample,
      cell,
      pre_qc_spots,
      post_qc_spots,
      nCount_Spatial,
      nFeature_Spatial,
      percent.mt,
      Glycolysis1,
      median_Glycolysis1,
      GlycoGroup,
      n_glyco_genes_present,
      glyco_genes_missing
    )

  spot_df <- meta_df %>%
    left_join(expr_df, by = "cell") %>%
    left_join(coords, by = "cell") %>%
    mutate(
      four_gene_score =
        four_gene_coef["TPI1"] * TPI1 +
        four_gene_coef["ENO1"] * ENO1 +
        four_gene_coef["LDHA"] * LDHA +
        four_gene_coef["SLC2A1"] * SLC2A1,
      MIF_SPP1_ligand_score = rowMeans(
        as.data.frame(select(., MIF, SPP1)),
        na.rm = TRUE
      )
    )

  spot_source_list[[sid]] <- spot_df

  sample_summary <- bind_rows(lapply(outcome_features, function(feat) {
    dat <- spot_df %>%
      filter(!is.na(.data[[feat]]), !is.na(GlycoGroup))

    high_vals <- dat %>%
      filter(GlycoGroup == "GlycoHigh") %>%
      pull(.data[[feat]])

    low_vals <- dat %>%
      filter(GlycoGroup == "GlycoLow") %>%
      pull(.data[[feat]])

    if (length(high_vals) == 0 || length(low_vals) == 0) {
      return(tibble(
        sample = sid,
        feature = feat,
        n_high = length(high_vals),
        n_low = length(low_vals),
        mean_high = NA_real_,
        mean_low = NA_real_,
        median_high = NA_real_,
        median_low = NA_real_,
        effect_high_minus_low = NA_real_,
        wilcox_p_greater = NA_real_
      ))
    }

    wt <- suppressWarnings(
      wilcox.test(high_vals, low_vals, alternative = "greater")
    )

    tibble(
      sample = sid,
      feature = feat,
      n_high = length(high_vals),
      n_low = length(low_vals),
      mean_high = mean(high_vals, na.rm = TRUE),
      mean_low = mean(low_vals, na.rm = TRUE),
      median_high = median(high_vals, na.rm = TRUE),
      median_low = median(low_vals, na.rm = TRUE),
      effect_high_minus_low = mean(high_vals, na.rm = TRUE) - mean(low_vals, na.rm = TRUE),
      wilcox_p_greater = wt$p.value
    )
  })) %>%
    group_by(sample) %>%
    mutate(wilcox_FDR_within_sample = p.adjust(wilcox_p_greater, method = "BH")) %>%
    ungroup()

  summary_list[[sid]] <- sample_summary

  group_counts <- table(seu$GlycoGroup)
  cat(
    sid,
    " GlycoHigh = ", as.integer(group_counts["GlycoHigh"]),
    "; GlycoLow = ", as.integer(group_counts["GlycoLow"]),
    "\n",
    sep = ""
  )

  # Candidate Figure 9A spatial feature plots for the representative sample.
  if (sid == representative_sample) {
    if (length(Images(seu)) > 0) {
      feature_plot_objects[["ENO1"]] <- SpatialFeaturePlot(
        seu,
        features = "ENO1",
        pt.size.factor = 1.4
      ) + ggtitle("ENO1")

      feature_plot_objects[["Glycolysis1"]] <- SpatialFeaturePlot(
        seu,
        features = "Glycolysis1",
        pt.size.factor = 1.4
      ) + ggtitle("Glycolysis1")

      feature_plot_objects[["SPP1"]] <- SpatialFeaturePlot(
        seu,
        features = "SPP1",
        pt.size.factor = 1.4
      ) + ggtitle("SPP1")

      feature_plot_objects[["MIF"]] <- SpatialFeaturePlot(
        seu,
        features = "MIF",
        pt.size.factor = 1.4
      ) + ggtitle("MIF")
    } else {
      feature_plot_objects[["ENO1"]] <- make_coordinate_feature_plot(spot_df, "ENO1", "ENO1")
      feature_plot_objects[["Glycolysis1"]] <- make_coordinate_feature_plot(spot_df, "Glycolysis1", "Glycolysis1")
      feature_plot_objects[["SPP1"]] <- make_coordinate_feature_plot(spot_df, "SPP1", "SPP1")
      feature_plot_objects[["MIF"]] <- make_coordinate_feature_plot(spot_df, "MIF", "MIF")
    }
  }

  # Supplementary Figure S18 group spatial plots.
  title_text <- sprintf(
    "%s (GlycoHigh = %d / GlycoLow = %d)",
    sid,
    as.integer(group_counts["GlycoHigh"]),
    as.integer(group_counts["GlycoLow"])
  )

  if (length(Images(seu)) > 0) {
    group_plot_objects[[sid]] <- SpatialDimPlot(
      seu,
      group.by = "GlycoGroup",
      cols = c("GlycoLow" = "#4575B4", "GlycoHigh" = "#D73027"),
      pt.size.factor = 1.4,
      alpha = c(0.85, 0.85)
    ) +
      ggtitle(title_text) +
      theme(
        plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
        legend.position = "right"
      )
  } else {
    group_plot_objects[[sid]] <- make_coordinate_group_plot(spot_df, title_text)
  }
}

spot_source <- bind_rows(spot_source_list)
spatial_summary <- bind_rows(summary_list) %>%
  group_by(feature) %>%
  mutate(wilcox_FDR_across_all_tests = p.adjust(wilcox_p_greater, method = "BH")) %>%
  ungroup()

source_file <- file.path(outdir, "Fig9_S18_GSE238264_spatial_spot_source_data.csv")
summary_file <- file.path(outdir, "Fig9_GSE238264_spatial_GlycoHigh_vs_GlycoLow_summary.csv")

readr::write_csv(spot_source, source_file)
readr::write_csv(spatial_summary, summary_file)

cat("Source CSV written to: ", source_file, "\n", sep = "")
cat("Summary CSV written to: ", summary_file, "\n", sep = "")

# -----------------------------
# 4. Optional RCTD-adjusted regression
# -----------------------------

rctd_results <- run_rctd_regression_if_available(
  source_df = spot_source,
  rctd_file = rctd_proportion_file,
  outdir = outdir
)

# -----------------------------
# 5. QC table
# -----------------------------

observed_counts <- spot_source %>%
  count(sample, GlycoGroup, name = "n") %>%
  tidyr::pivot_wider(names_from = GlycoGroup, values_from = n, values_fill = 0) %>%
  rename(
    observed_low = GlycoLow,
    observed_high = GlycoHigh
  ) %>%
  left_join(expected_group_counts, by = "sample") %>%
  mutate(
    exact_group_count_pass =
      observed_high == expected_high & observed_low == expected_low
  )

feature_presence_by_sample <- spot_source %>%
  group_by(sample) %>%
  summarise(
    missing_feature_columns = paste(
      outcome_features[!outcome_features %in% names(cur_data())],
      collapse = ";"
    ),
    .groups = "drop"
  )

direction_summary <- spatial_summary %>%
  group_by(feature) %>%
  summarise(
    n_samples = n_distinct(sample),
    n_positive = sum(effect_high_minus_low > 0, na.rm = TRUE),
    all_positive = all(effect_high_minus_low > 0, na.rm = TRUE),
    median_effect = median(effect_high_minus_low, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(expected_median_effects, by = "feature") %>%
  mutate(
    reference_tolerance_pass = ifelse(
      is.na(expected_median_effect),
      TRUE,
      abs(median_effect - expected_median_effect) <= 0.20
    )
  )

qc <- bind_rows(
  qc_row(
    "all_samples_processed",
    paste(sort(unique(spot_source$sample)), collapse = "; "),
    paste(sample_ids, collapse = "; "),
    setequal(unique(spot_source$sample), sample_ids)
  ),
  qc_row(
    "total_primary_spatial_spots_after_QC",
    nrow(spot_source),
    sum(expected_group_counts$expected_high + expected_group_counts$expected_low),
    nrow(spot_source) == sum(expected_group_counts$expected_high + expected_group_counts$expected_low)
  ),
  bind_rows(lapply(seq_len(nrow(observed_counts)), function(i) {
    qc_row(
      paste0("group_counts_", observed_counts$sample[i]),
      paste0(
        "GlycoHigh=", observed_counts$observed_high[i],
        "; GlycoLow=", observed_counts$observed_low[i]
      ),
      paste0(
        "GlycoHigh=", observed_counts$expected_high[i],
        "; GlycoLow=", observed_counts$expected_low[i]
      ),
      observed_counts$exact_group_count_pass[i]
    )
  })),
  qc_row(
    "all_required_outcome_features_available",
    paste(outcome_features[outcome_features %in% names(spot_source)], collapse = "; "),
    paste(outcome_features, collapse = "; "),
    all(outcome_features %in% names(spot_source))
  ),
  bind_rows(lapply(outcome_features, function(feat) {
    ds <- direction_summary %>% filter(feature == feat)

    qc_row(
      paste0("direction_positive_all_samples_", feat),
      ifelse(
        nrow(ds) == 0,
        "missing",
        paste0("n_positive=", ds$n_positive, "/", ds$n_samples, "; median_effect=", signif(ds$median_effect, 4))
      ),
      "4/4 samples positive",
      nrow(ds) == 1 && ds$n_samples == 4 && ds$n_positive == 4
    )
  })),
  bind_rows(lapply(seq_len(nrow(direction_summary)), function(i) {
    qc_row(
      paste0("reference_median_effect_", direction_summary$feature[i]),
      signif(direction_summary$median_effect[i], 4),
      direction_summary$expected_median_effect[i],
      direction_summary$reference_tolerance_pass[i],
      severity = "CHECK",
      notes = "Reference tolerance check only. Mandatory output gating uses group counts and 4/4 positive direction consistency."
    )
  })),
  qc_row(
    "optional_RCTD_file_present",
    file.exists(rctd_proportion_file),
    paste0("Optional: ", rctd_proportion_file),
    TRUE,
    severity = "INFO",
    notes = "If absent, primary spatial source/QC/figures still run. RCTD-adjusted Figure S20 can be regenerated by the dedicated restoration script."
  ),
  qc_row(
    "optional_RCTD_results_written",
    !is.null(rctd_results),
    "TRUE if rctd_spot_proportions.csv was supplied",
    TRUE,
    severity = "INFO",
    notes = "Informational only."
  )
)

qc_file <- file.path(outdir, "Fig9_S18_GSE238264_spatial_QC_check.csv")
readr::write_csv(qc, qc_file)

cat("QC CSV written to: ", qc_file, "\n", sep = "")
print(qc)

mandatory_pass <- all(qc$pass[qc$severity == "ERROR"])

if (!mandatory_pass) {
  stop(
    "Mandatory QC failed. Source, summary, and QC CSV files were written, ",
    "but PDF/PNG files were not generated. Inspect: ",
    qc_file
  )
}

cat("Mandatory QC PASS. Proceeding to PDF/PNG generation.\n")

# -----------------------------
# 6. Figure 9A: representative spatial feature maps
# -----------------------------

if (length(feature_plot_objects) == 4) {
  fig9a <- (
    feature_plot_objects[["ENO1"]] |
      feature_plot_objects[["Glycolysis1"]]
  ) / (
    feature_plot_objects[["SPP1"]] |
      feature_plot_objects[["MIF"]]
  )

  fig9a_pdf <- file.path(outdir, paste0("Fig9A_spatial_feature_maps_", representative_sample, ".pdf"))
  fig9a_png <- file.path(outdir, paste0("Fig9A_spatial_feature_maps_", representative_sample, ".png"))

  ggsave(
    filename = fig9a_pdf,
    plot = fig9a,
    width = 10,
    height = 8,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = fig9a_png,
    plot = fig9a,
    width = 10,
    height = 8,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  cat("Figure 9A PDF written to: ", fig9a_pdf, "\n", sep = "")
  cat("Figure 9A PNG written to: ", fig9a_png, "\n", sep = "")
} else {
  warning("Figure 9A feature plot objects were incomplete; Fig9A was not generated.")
}

# -----------------------------
# 7. Figure 9B: spatial direction-consistency summary
# -----------------------------

feature_labels <- c(
  ENO1 = "ENO1",
  four_gene_score = "Four-gene score",
  PTGES = "PTGES",
  MIF = "MIF",
  MIF_SPP1_ligand_score = "MIF/SPP1 score",
  SPP1 = "SPP1"
)

plot_summary <- spatial_summary %>%
  mutate(
    feature = factor(feature, levels = outcome_features),
    feature_label = feature_labels[as.character(feature)],
    direction = ifelse(effect_high_minus_low > 0, "Higher in GlycoHigh", "Not higher in GlycoHigh"),
    neg_log10_p = -log10(pmax(wilcox_p_greater, .Machine$double.xmin)),
    sig_label = case_when(
      wilcox_FDR_across_all_tests < 0.001 ~ "***",
      wilcox_FDR_across_all_tests < 0.01 ~ "**",
      wilcox_FDR_across_all_tests < 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

fig9b <- ggplot(
  plot_summary,
  aes(
    x = sample,
    y = feature_label,
    size = neg_log10_p,
    color = effect_high_minus_low
  )
) +
  geom_point(alpha = 0.9) +
  geom_text(aes(label = sig_label), vjust = -1.1, size = 3.2) +
  scale_size_continuous(name = "-log10(p)", range = c(2.5, 11)) +
  scale_color_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#D73027",
    midpoint = 0,
    name = "High - Low"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 11),
    legend.position = "right",
    plot.title = element_blank()
  )

fig9b_pdf <- file.path(outdir, "Fig9B_spatial_GlycoHigh_vs_GlycoLow_direction_consistency.pdf")
fig9b_png <- file.path(outdir, "Fig9B_spatial_GlycoHigh_vs_GlycoLow_direction_consistency.png")

ggsave(
  filename = fig9b_pdf,
  plot = fig9b,
  width = 8.5,
  height = 5.2,
  units = "in",
  bg = "white"
)

ggsave(
  filename = fig9b_png,
  plot = fig9b,
  width = 8.5,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

cat("Figure 9B PDF written to: ", fig9b_pdf, "\n", sep = "")
cat("Figure 9B PNG written to: ", fig9b_png, "\n", sep = "")

# -----------------------------
# 8. Supplementary Figure S18: spatial GlycoHigh/GlycoLow maps
# -----------------------------

if (length(group_plot_objects) == length(sample_ids)) {
  figs18 <- patchwork::wrap_plots(group_plot_objects, ncol = 2)

  figs18_pdf <- file.path(outdir, "FigS18_spatial_GlycoHigh_GlycoLow_spot_stratification.pdf")
  figs18_png <- file.path(outdir, "FigS18_spatial_GlycoHigh_GlycoLow_spot_stratification.png")

  ggsave(
    filename = figs18_pdf,
    plot = figs18,
    width = 12,
    height = 10,
    units = "in",
    bg = "white"
  )

  ggsave(
    filename = figs18_png,
    plot = figs18,
    width = 12,
    height = 10,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  cat("Supplementary Figure S18 PDF written to: ", figs18_pdf, "\n", sep = "")
  cat("Supplementary Figure S18 PNG written to: ", figs18_png, "\n", sep = "")
} else {
  warning("S18 group plot objects were incomplete; FigS18 was not generated.")
}

# -----------------------------
# 9. Console summary
# -----------------------------

cat("\nSpatial transcriptomics analysis complete.\n")
cat("Interpretation reminder: Visium results are spot-level tissue co-enrichment analyses only.\n")
cat("Source CSV: ", source_file, "\n", sep = "")
cat("Summary CSV: ", summary_file, "\n", sep = "")
cat("QC CSV: ", qc_file, "\n", sep = "")

cat("\nObserved GlycoHigh/GlycoLow counts:\n")
print(observed_counts)

cat("\nDirection consistency summary:\n")
print(direction_summary)
