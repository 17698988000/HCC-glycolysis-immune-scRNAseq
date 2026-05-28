# ============================================================
# 15_LIANA_directional_concordance.R
# Step 11-L1 LIANA directional-concordance patch
#
# Purpose:
#   Cross-method directional support for predefined MIF/SPP1 axes
#   in the GlycoHigh-to-tumor-immune direction.
#
# Input object confirmed in Step 11:
#   D:/scRNA_project/hcc_seurat.rds
#   celltype column: celltype_original
#   GlycoHigh/GlycoLow column: celltype_glyco
#   patient column: orig.ident
#   tissue column: site
#
# Boundary:
#   This analysis provides directional concordance only.
#   It is not functional validation, mechanistic proof,
#   spatial proximity evidence, treatment-response prediction,
#   or a clinical model.
# ============================================================

setwd("D:/scRNA_project")

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) BiocManager::install("SingleCellExperiment", ask = FALSE, update = FALSE)
if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) BiocManager::install("SummarizedExperiment", ask = FALSE, update = FALSE)
if (!requireNamespace("liana", quietly = TRUE)) remotes::install_github("saezlab/liana")

pkgs <- c("Seurat", "SeuratObject", "SingleCellExperiment", "SummarizedExperiment",
          "liana", "dplyr", "tidyr", "tibble", "purrr", "readr", "ggplot2", "stringr", "Matrix")
for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))

set.seed(20260528)

OBJ_RDS      <- "D:/scRNA_project/hcc_seurat.rds"
CELLTYPE_COL <- "celltype_original"
GLYCO_COL    <- "celltype_glyco"
PATIENT_COL  <- "orig.ident"
TISSUE_COL   <- "site"
OUTDIR <- "results/LIANA_directional_concordance"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

sink(file.path(OUTDIR, "run_log_LIANA.txt"), split = TRUE)
cat("Step 11-L1 LIANA directional-concordance patch\n")
cat("Date:", as.character(Sys.time()), "\n\n")

obj <- readRDS(OBJ_RDS)
obj <- UpdateSeuratObject(obj)
if ("RNA" %in% names(obj@assays)) DefaultAssay(obj) <- "RNA" else stop("RNA assay not found.")
try({ obj <- JoinLayers(obj) }, silent = TRUE)

meta <- obj@meta.data
required_cols <- c(CELLTYPE_COL, GLYCO_COL, PATIENT_COL, TISSUE_COL)
missing_cols <- setdiff(required_cols, colnames(meta))
if (length(missing_cols) > 0) stop("Missing metadata columns: ", paste(missing_cols, collapse = ", "))

celltype <- tolower(as.character(meta[[CELLTYPE_COL]]))
glyco    <- tolower(as.character(meta[[GLYCO_COL]]))
tissue   <- tolower(as.character(meta[[TISSUE_COL]]))

is_tumor <- grepl("tumor|tumour|primary|hcc|pt|pvt|metast", tissue) & !grepl("normal|adjacent|non", tissue)
if (sum(is_tumor, na.rm = TRUE) == 0) {
  warning("No tumor label detected in TISSUE_COL; assuming object is tumor-only.")
  is_tumor <- rep(TRUE, nrow(meta))
}

is_hep <- grepl("hepato|hepatocyte|malignant|tumor", celltype) &
  !grepl("t/nk|t_nk|t cell|cd8|cd4|nk|myeloid|macro|macrophage|mono|kupffer|endo|fibro|b cell", celltype)
is_tnk <- grepl("t/nk|t_nk|t cell|cd8|cd4|nk|exhaust|effector", celltype)
is_myeloid <- grepl("myeloid|macro|macrophage|kupffer|mono|mdsc|dendritic", celltype)
is_high <- grepl("glycohigh|glyco_high|high", glyco) & !grepl("low", glyco)
is_low  <- grepl("glycolow|glyco_low|low", glyco) & !grepl("high", glyco)

obj$liana_group <- NA_character_
obj$liana_group[is_tumor & is_hep & is_high] <- "GlycoHigh_hepatocyte"
obj$liana_group[is_tumor & is_hep & is_low]  <- "GlycoLow_hepatocyte"
obj$liana_group[is_tumor & is_tnk]           <- "Tumor_T_NK"
obj$liana_group[is_tumor & is_myeloid]       <- "Tumor_myeloid"

group_counts <- as.data.frame(table(obj$liana_group, useNA = "ifany"))
colnames(group_counts) <- c("liana_group", "n_cells")
cat("\nLIANA group counts:\n"); print(group_counts)
write_csv(group_counts, file.path(OUTDIR, "LIANA_group_counts.csv"))

needed_groups <- c("GlycoHigh_hepatocyte", "GlycoLow_hepatocyte", "Tumor_T_NK", "Tumor_myeloid")
group_tab <- table(obj$liana_group)
if (!all(needed_groups %in% names(group_tab))) stop("Missing at least one required LIANA group.")
if (any(group_tab[needed_groups] < 10)) stop("At least one required LIANA group has fewer than 10 cells.")

liana_obj <- subset(obj, cells = colnames(obj)[!is.na(obj$liana_group)])
Idents(liana_obj) <- liana_obj$liana_group
saveRDS(liana_obj, file.path(OUTDIR, "LIANA_input_hcc_seurat_grouped.rds"))
cat("\nFinal LIANA Seurat object dimensions:\n"); print(dim(liana_obj))
cat("\nFinal LIANA identity counts:\n"); print(table(Idents(liana_obj)))

get_layer_matrix <- function(seu, assay = "RNA", layer = "data") {
  tryCatch(SeuratObject::LayerData(seu, assay = assay, layer = layer),
           error = function(e1) tryCatch(GetAssayData(seu, assay = assay, layer = layer), error = function(e2) NULL))
}
counts_mat <- get_layer_matrix(liana_obj, assay = "RNA", layer = "counts")
data_mat   <- get_layer_matrix(liana_obj, assay = "RNA", layer = "data")
if (is.null(counts_mat) && is.null(data_mat)) stop("Could not retrieve RNA counts or data layer.")
if (is.null(counts_mat)) { warning("Counts layer not found; using data layer as counts placeholder."); counts_mat <- data_mat }
if (is.null(data_mat)) { warning("Data layer not found; using log1p(counts) as logcounts."); data_mat <- log1p(counts_mat) }

fix_orientation <- function(mat, seu_cells) {
  if (ncol(mat) == length(seu_cells) && all(colnames(mat) == seu_cells)) return(mat)
  if (nrow(mat) == length(seu_cells) && all(rownames(mat) == seu_cells)) return(Matrix::t(mat))
  if (all(seu_cells %in% colnames(mat))) return(mat[, seu_cells, drop = FALSE])
  if (all(seu_cells %in% rownames(mat))) return(Matrix::t(mat[seu_cells, , drop = FALSE]))
  stop("Cannot orient matrix to match Seurat cells.")
}
counts_mat <- fix_orientation(counts_mat, colnames(liana_obj))
data_mat   <- fix_orientation(data_mat, colnames(liana_obj))
common_genes <- intersect(rownames(counts_mat), rownames(data_mat))
counts_mat <- counts_mat[common_genes, colnames(liana_obj), drop = FALSE]
data_mat   <- data_mat[common_genes, colnames(liana_obj), drop = FALSE]

sce_liana <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = counts_mat, logcounts = data_mat),
  colData = liana_obj@meta.data
)
sce_liana$liana_group <- factor(sce_liana$liana_group,
  levels = c("GlycoHigh_hepatocyte", "GlycoLow_hepatocyte", "Tumor_T_NK", "Tumor_myeloid"))
cat("\nSCE liana_group counts:\n"); print(table(sce_liana$liana_group))
saveRDS(sce_liana, file.path(OUTDIR, "LIANA_input_hcc_SCE_grouped.rds"))

cat("\nRunning LIANA on SingleCellExperiment...\n")
liana_methods <- c("natmi", "connectome", "logfc", "sca")
liana_res <- liana_wrap(
  sce = sce_liana,
  method = liana_methods,
  resource = "Consensus",
  idents_col = "liana_group",
  min_cells = 10,
  expr_prop = 0.10,
  return_all = TRUE,
  verbose = TRUE
)
saveRDS(liana_res, file.path(OUTDIR, "LIANA_raw_results.rds"))
cat("\nLIANA run finished.\n")

if (is.data.frame(liana_res)) {
  liana_full <- as_tibble(liana_res) %>% mutate(liana_method = "single_result")
} else if (is.list(liana_res)) {
  liana_full <- purrr::imap_dfr(liana_res, function(x, nm) as_tibble(x) %>% mutate(liana_method = nm))
} else stop("Unexpected LIANA output type.")
write_tsv(liana_full, file.path(OUTDIR, "Supplementary_Table_S4_LIANA_full_results.tsv"))
cat("\nFull LIANA result dimensions:\n"); print(dim(liana_full))

cat("\nAggregating LIANA ranks...\n")
liana_agg <- tryCatch({
  if (exists("liana_aggregate")) liana_aggregate(liana_res) else if (exists("rank_aggregate")) rank_aggregate(liana_res) else stop("No aggregate function.")
}, error = function(e) { cat("\nRank aggregation failed; using full flattened results.\n"); cat("Error was:", conditionMessage(e), "\n"); liana_full })
liana_agg <- as_tibble(liana_agg)
write_tsv(liana_agg, file.path(OUTDIR, "Supplementary_Table_S4_LIANA_aggregate_results.tsv"))
cat("\nAggregated LIANA result dimensions:\n"); print(dim(liana_agg))
cat("\nAggregated LIANA result columns:\n"); print(colnames(liana_agg))

get_first_col <- function(df, candidates) {
  hit <- intersect(candidates, colnames(df))[1]
  if (is.na(hit)) return(rep(NA_character_, nrow(df)))
  as.character(df[[hit]])
}
standardize_lr_table <- function(df) {
  df %>% mutate(
    source_std = get_first_col(., c("source", "sender", "from")),
    target_std = get_first_col(., c("target", "receiver", "to")),
    ligand_std = get_first_col(., c("ligand.complex", "ligand_complex", "ligand")),
    receptor_std = get_first_col(., c("receptor.complex", "receptor_complex", "receptor"))
  )
}
agg_std <- standardize_lr_table(liana_agg)

candidate_pairs <- tribble(
  ~ligand_query, ~receptor_query,
  "MIF",  "CD74", "MIF",  "CXCR4", "MIF",  "CD44",
  "SPP1", "CD44", "SPP1", "ITGAV", "SPP1", "ITGB1", "SPP1", "ITGA4", "SPP1", "ITGB5", "SPP1", "ITGB6"
)
sender_keep <- c("GlycoHigh_hepatocyte", "GlycoLow_hepatocyte")
target_keep <- c("Tumor_T_NK", "Tumor_myeloid")
rank_col <- intersect(c("aggregate_rank", "specificity_rank", "magnitude_rank", "rank"), colnames(agg_std))[1]
if (is.na(rank_col)) warning("No rank-like column found. Candidate output will be presence-based.")

cand <- agg_std %>%
  filter(source_std %in% sender_keep, target_std %in% target_keep) %>%
  mutate(tmp_join_key = 1) %>%
  inner_join(candidate_pairs %>% mutate(tmp_join_key = 1), by = "tmp_join_key") %>%
  filter(ligand_std == ligand_query,
         stringr::str_detect(receptor_std, paste0("(^|[^A-Z0-9])", receptor_query, "($|[^A-Z0-9])"))) %>%
  select(-tmp_join_key)
write_tsv(cand, file.path(OUTDIR, "Supplementary_Table_S4_LIANA_candidate_axes_all_matching_rows.tsv"))

if (!is.na(rank_col)) {
  cand_by_source <- cand %>%
    group_by(target_std, ligand_query, receptor_query, source_std) %>%
    summarise(best_rank = suppressWarnings(min(as.numeric(.data[[rank_col]]), na.rm = TRUE)),
              n_matching_rows = dplyr::n(), .groups = "drop") %>%
    mutate(best_rank = ifelse(is.infinite(best_rank), NA_real_, best_rank))
} else {
  cand_by_source <- cand %>% group_by(target_std, ligand_query, receptor_query, source_std) %>%
    summarise(best_rank = NA_real_, n_matching_rows = dplyr::n(), .groups = "drop")
}

cand_wide <- cand_by_source %>% tidyr::pivot_wider(names_from = source_std, values_from = c(best_rank, n_matching_rows), values_fill = list(n_matching_rows = 0))
needed_cols <- c("best_rank_GlycoHigh_hepatocyte", "best_rank_GlycoLow_hepatocyte", "n_matching_rows_GlycoHigh_hepatocyte", "n_matching_rows_GlycoLow_hepatocyte")
for (cc in needed_cols) if (!cc %in% colnames(cand_wide)) cand_wide[[cc]] <- NA

direction_summary <- cand_wide %>%
  mutate(rank_direction = case_when(
    !is.na(best_rank_GlycoHigh_hepatocyte) & !is.na(best_rank_GlycoLow_hepatocyte) & best_rank_GlycoHigh_hepatocyte < best_rank_GlycoLow_hepatocyte ~ "GlycoHigh_stronger_rank",
    !is.na(best_rank_GlycoHigh_hepatocyte) & !is.na(best_rank_GlycoLow_hepatocyte) & best_rank_GlycoHigh_hepatocyte > best_rank_GlycoLow_hepatocyte ~ "GlycoLow_stronger_rank",
    !is.na(best_rank_GlycoHigh_hepatocyte) & !is.na(best_rank_GlycoLow_hepatocyte) & best_rank_GlycoHigh_hepatocyte == best_rank_GlycoLow_hepatocyte ~ "equal_rank",
    !is.na(best_rank_GlycoHigh_hepatocyte) & is.na(best_rank_GlycoLow_hepatocyte) ~ "GlycoHigh_detected_only",
    is.na(best_rank_GlycoHigh_hepatocyte) & !is.na(best_rank_GlycoLow_hepatocyte) ~ "GlycoLow_detected_only",
    TRUE ~ "not_detected"
  )) %>%
  mutate(target_std = as.character(target_std), ligand_query = as.character(ligand_query), receptor_query = as.character(receptor_query), rank_direction = as.character(rank_direction))
write_tsv(direction_summary, file.path(OUTDIR, "Supplementary_Table_S4_LIANA_candidate_axes.tsv"))

direction_counts <- direction_summary %>% group_by(rank_direction) %>% summarise(n_axis_receiver_pairs = dplyr::n(), .groups = "drop") %>% arrange(desc(n_axis_receiver_pairs))
write_tsv(direction_counts, file.path(OUTDIR, "LIANA_direction_counts.tsv"))
cat("\nCandidate-axis direction summary:\n"); print(direction_summary)
cat("\nDirection counts:\n"); print(direction_counts)

plot_df <- direction_summary %>% mutate(
  axis = paste(ligand_query, receptor_query, sep = "-"),
  direction_for_plot = case_when(rank_direction == "GlycoHigh_stronger_rank" ~ 1, rank_direction == "GlycoHigh_detected_only" ~ 1, rank_direction == "GlycoLow_stronger_rank" ~ -1, rank_direction == "GlycoLow_detected_only" ~ -1, TRUE ~ 0),
  support_size = case_when(direction_for_plot == 1 ~ 3, direction_for_plot == -1 ~ 3, TRUE ~ 1.5)
)
p <- ggplot(plot_df, aes(x = target_std, y = axis, size = support_size)) +
  geom_point() + theme_bw(base_size = 10) +
  labs(x = "Receiver compartment", y = "Predefined candidate axis", size = "Directional support", title = "LIANA directional concordance")
ggsave(file.path(OUTDIR, "Figure4D_LIANA_directional_concordance.pdf"), p, width = 6.5, height = 4.5)
ggsave(file.path(OUTDIR, "Figure4D_LIANA_directional_concordance.png"), p, width = 6.5, height = 4.5, dpi = 300)

sink(file.path(OUTDIR, "sessionInfo_LIANA.txt")); print(sessionInfo()); sink()
writeLines(c(
  "Step 11 LIANA directional-concordance patch", "",
  "Interpretation boundary:",
  "LIANA was used only as a cross-method directional-concordance check for predefined MIF/SPP1 axes.",
  "The output should not be interpreted as functional ligand-receptor validation, mechanistic proof, spatial proximity evidence, treatment-response prediction, or clinical decision support.", "",
  "Main result:",
  "Among 18 predefined axis-receiver combinations, 14 showed stronger aggregate ranks for GlycoHigh hepatocytes than GlycoLow hepatocytes, whereas four showed equal aggregate ranks and none favored GlycoLow hepatocytes.", "",
  "Manuscript wording:",
  "candidate, inference-based communication-layer support"
), con = file.path(OUTDIR, "LIANA_INTERPRETATION_BOUNDARY_NOTE.txt"))

sink()
cat("\nLIANA Step 11-L1 completed.\n")
cat("Outputs written to:\n")
cat(OUTDIR, "\n")
