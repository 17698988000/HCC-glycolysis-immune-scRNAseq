#!/usr/bin/env Rscript

# Spatial-block permutation inference for composition-adjusted Visium models.
#
# Required inputs:
#   SPATIAL_SOURCE_CSV: Fig9_S18_GSE238264_spatial_spot_source_data.csv
#   RCTD_PROPORTION_CSV: rctd_spot_proportions.csv
#
# The script fits section-specific composition-adjusted models and estimates
# empirical two-sided P values by shuffling the glycolysis predictor within
# spatial grid blocks. This preserves broad section-level spatial structure
# more closely than unrestricted spot permutation.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
})

set.seed(20260603)

source_csv <- Sys.getenv(
  "SPATIAL_SOURCE_CSV",
  unset = file.path("results", "GSE238264_spatial", "Fig9_S18_GSE238264_spatial_spot_source_data.csv")
)
rctd_csv <- Sys.getenv(
  "RCTD_PROPORTION_CSV",
  unset = file.path("data", "GSE238264", "rctd_spot_proportions.csv")
)
outdir <- Sys.getenv(
  "SPATIAL_PERMUTATION_OUTDIR",
  unset = file.path("results", "GSE238264_spatial_block_permutation")
)
n_permutations <- as.integer(Sys.getenv("N_SPATIAL_PERMUTATIONS", unset = "999"))
n_x_blocks <- as.integer(Sys.getenv("N_X_BLOCKS", unset = "4"))
n_y_blocks <- as.integer(Sys.getenv("N_Y_BLOCKS", unset = "4"))
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_csv)) stop("Missing spatial source CSV: ", source_csv, call. = FALSE)
if (!file.exists(rctd_csv)) stop("Missing RCTD proportion CSV: ", rctd_csv, call. = FALSE)

outcomes <- c("ENO1", "four_gene_score", "PTGES", "MIF", "SPP1", "MIF_SPP1_ligand_score")

standardize_rctd_cols <- function(df) {
  names(df) <- tolower(str_replace_all(names(df), "[^a-z0-9]+", "_"))
  sample_col <- intersect(c("sample", "sample_id", "sid"), names(df))[1]
  cell_col <- intersect(c("cell", "barcode", "spot", "spot_barcode"), names(df))[1]
  if (is.na(sample_col) || is.na(cell_col)) stop("RCTD file needs sample and cell/barcode columns.")
  df <- df %>% rename(sample = all_of(sample_col), cell = all_of(cell_col))
  copy_one <- function(candidates, target) {
    hit <- intersect(candidates, names(df))[1]
    if (is.na(hit)) stop("Missing RCTD proportion column for ", target, call. = FALSE)
    df[[target]] <<- as.numeric(df[[hit]])
  }
  copy_one(c("hepatocyte", "hepatocytes"), "prop_hepatocyte")
  copy_one(c("myeloid", "myeloid_cells", "macrophage", "macrophages"), "prop_myeloid")
  copy_one(c("t_nk", "tnk", "t_nk_cells", "t_cells_nk_cells"), "prop_t_nk")
  copy_one(c("fibroblast", "fibroblasts"), "prop_fibroblast")
  copy_one(c("endothelial", "endothelial_cells"), "prop_endothelial")
  df %>% select(
    sample, cell, prop_hepatocyte, prop_myeloid, prop_t_nk,
    prop_fibroblast, prop_endothelial
  )
}

choose_coordinates <- function(df) {
  candidates <- list(
    c("imagecol", "imagerow"),
    c("array_col", "array_row"),
    c("x", "y")
  )
  for (pair in candidates) {
    if (all(pair %in% names(df))) return(pair)
  }
  stop("Spatial source CSV must contain imagecol/imagerow, array_col/array_row, or x/y coordinates.")
}

fit_beta <- function(dat, outcome, predictor = "Glycolysis1_z") {
  dat$outcome_z <- as.numeric(scale(dat[[outcome]]))
  formula <- as.formula(
    paste(
      "outcome_z ~", predictor,
      "+ prop_hepatocyte + prop_myeloid + prop_t_nk +",
      "prop_fibroblast + prop_endothelial + log_depth"
    )
  )
  fit <- lm(formula, data = dat)
  unname(coef(fit)[predictor])
}

source_df <- read_csv(source_csv, show_col_types = FALSE)
rctd <- read_csv(rctd_csv, show_col_types = FALSE) %>% standardize_rctd_cols()
coord_cols <- choose_coordinates(source_df)

merged <- source_df %>%
  left_join(rctd, by = c("sample", "cell")) %>%
  filter(
    !is.na(Glycolysis1),
    !is.na(nCount_Spatial),
    if_all(
      c(prop_hepatocyte, prop_myeloid, prop_t_nk, prop_fibroblast, prop_endothelial),
      ~ !is.na(.x)
    )
  ) %>%
  group_by(sample) %>%
  mutate(
    Glycolysis1_z = as.numeric(scale(Glycolysis1)),
    log_depth = log1p(nCount_Spatial),
    x_block = ntile(.data[[coord_cols[1]]], n_x_blocks),
    y_block = ntile(.data[[coord_cols[2]]], n_y_blocks),
    spatial_block = interaction(x_block, y_block, drop = TRUE)
  ) %>%
  ungroup()

results <- list()
null_results <- list()
for (sid in sort(unique(merged$sample))) {
  for (outcome in outcomes) {
    dat <- merged %>% filter(sample == sid, !is.na(.data[[outcome]]))
    observed_beta <- fit_beta(dat, outcome)
    null_beta <- replicate(n_permutations, {
      permuted <- dat %>%
        group_by(spatial_block) %>%
        mutate(Glycolysis1_z_perm = sample(Glycolysis1_z, replace = FALSE)) %>%
        ungroup()
      fit_beta(permuted, outcome, "Glycolysis1_z_perm")
    })
    empirical_p <- (1 + sum(abs(null_beta) >= abs(observed_beta), na.rm = TRUE)) /
      (1 + sum(!is.na(null_beta)))
    results[[length(results) + 1]] <- tibble(
      sample = sid,
      outcome = outcome,
      n_spots = nrow(dat),
      n_spatial_blocks = n_distinct(dat$spatial_block),
      observed_beta = observed_beta,
      empirical_p_two_sided = empirical_p,
      n_permutations = n_permutations
    )
    null_results[[length(null_results) + 1]] <- tibble(
      sample = sid,
      outcome = outcome,
      permutation = seq_along(null_beta),
      null_beta = null_beta
    )
  }
}

result_df <- bind_rows(results) %>%
  group_by(outcome) %>%
  mutate(empirical_FDR_within_outcome = p.adjust(empirical_p_two_sided, method = "BH")) %>%
  ungroup() %>%
  mutate(empirical_FDR_global = p.adjust(empirical_p_two_sided, method = "BH"))
null_df <- bind_rows(null_results)

write_csv(result_df, file.path(outdir, "spatial_block_permutation_section_results.csv"))
write_csv(null_df, file.path(outdir, "spatial_block_permutation_null_betas.csv"))
message("Spatial-block permutation outputs written to: ", outdir)
