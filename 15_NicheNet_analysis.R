# =============================================================================
# 15_NicheNet_analysis.R
# NicheNet Ligand Activity Analysis (Section 2.16)
# Prioritizes ligands from GlycoHigh hepatocytes predicted to regulate
# immunosuppression and T cell exhaustion markers in T/NK and myeloid receivers
# Generates: Supplementary Figure S22
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(nichenetr)   # v2.2.1
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(RColorBrewer)
})

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Load processed Seurat object and NicheNet reference databases
# -----------------------------------------------------------------------------
message("Loading data...")
seurat_obj <- readRDS("results/seurat_hcc_annotated.rds")

# NicheNet databases from Zenodo (https://zenodo.org/record/7074291)
# 脚本会自动检查并下载，无需手动操作

db_dir <- "data/nichenet_databases"
if (!dir.exists(db_dir)) dir.create(db_dir, recursive = TRUE)

# 数据库文件列表
db_files <- list(
  ligand_target_matrix = "ligand_target_matrix_nsga2r_final.rds",
  lr_network           = "lr_network_human_21122021.rds",
  weighted_networks    = "weighted_networks_nsga2r_final.rds"
)
zenodo_base <- "https://zenodo.org/record/7074291/files"

# 自动下载缺失的数据库文件
for (key in names(db_files)) {
  fpath <- file.path(db_dir, db_files[[key]])
  if (!file.exists(fpath)) {
    url <- paste0(zenodo_base, "/", db_files[[key]])
    message(sprintf("Downloading %s ...", db_files[[key]]))
    tryCatch(
      download.file(url, destfile = fpath, mode = "wb", method = "auto"),
      error = function(e) stop(sprintf(
        "下载失败: %s\n请手动从以下地址下载后放入 %s 目录:\n%s",
        db_files[[key]], db_dir, url
      ))
    )
    message(sprintf("  已保存至 %s", fpath))
  } else {
    message(sprintf("  %s 已存在，跳过下载", db_files[[key]]))
  }
}

# 加载数据库
message("Loading NicheNet databases...")
ligand_target_matrix <- readRDS(file.path(db_dir, db_files$ligand_target_matrix))
lr_network           <- readRDS(file.path(db_dir, db_files$lr_network))
weighted_networks    <- readRDS(file.path(db_dir, db_files$weighted_networks))
message("NicheNet databases loaded successfully.")

# -----------------------------------------------------------------------------
# 2. Define sender and receiver populations (Section 2.16)
# -----------------------------------------------------------------------------
# Sender: GlycoHigh hepatocytes (tumor-derived hepatocytes above median AUCell score)
# n = 5,589 per manuscript (subset of 15,391 tumor hepatocytes)

# Glycolysis AUC must be pre-computed (from Script 02)
if (!"Glycolysis_AUC" %in% colnames(seurat_obj@meta.data)) {
  stop("Glycolysis_AUC not found in metadata. Run 02_glycolysis_scoring.R first.")
}

tumor_hep <- subset(seurat_obj, subset = cell_type == "Hepatocyte" & tissue == "Tumor")
glyco_median <- median(tumor_hep@meta.data$Glycolysis_AUC, na.rm = TRUE)
tumor_hep$glyco_group <- ifelse(tumor_hep@meta.data$Glycolysis_AUC >= glyco_median,
                                 "GlycoHigh", "GlycoLow")

glyco_high_cells <- Cells(tumor_hep)[tumor_hep$glyco_group == "GlycoHigh"]
message(sprintf("GlycoHigh sender cells: %d (manuscript: 5,589)", length(glyco_high_cells)))

# Receiver: T/NK and myeloid cells (tumor-derived)
receiver_cells <- Cells(subset(seurat_obj,
  subset = cell_type %in% c("T_NK", "Myeloid", "CD8T") & tissue == "Tumor"))
message(sprintf("Receiver cells (T/NK + Myeloid): %d (manuscript: 11,383)", length(receiver_cells)))

# Extract expression matrices
sender_obj   <- subset(seurat_obj, cells = glyco_high_cells)
receiver_obj <- subset(seurat_obj, cells = receiver_cells)

# -----------------------------------------------------------------------------
# 3. Define expressed genes in each population
# Gene expressed = detected in >= 10% of cells within that population
# -----------------------------------------------------------------------------
message("Identifying expressed genes...")

get_expressed_genes <- function(seurat_subset, pct_threshold = 0.10) {
  expr_mat  <- GetAssayData(seurat_subset, slot = "counts")
  pct_expr  <- rowMeans(expr_mat > 0)
  names(pct_expr)[pct_expr >= pct_threshold]
}

expressed_genes_sender   <- get_expressed_genes(sender_obj)
expressed_genes_receiver <- get_expressed_genes(receiver_obj)

message(sprintf("Expressed genes — Sender: %d | Receiver: %d",
                length(expressed_genes_sender), length(expressed_genes_receiver)))

# -----------------------------------------------------------------------------
# 4. Define gene set of interest (immunosuppression / T cell exhaustion markers)
# 20 markers as in Section 2.16; subset to those expressed in receiver cells
# -----------------------------------------------------------------------------
geneset_of_interest_all <- c(
  "PDCD1", "HAVCR2", "TIGIT", "LAG3", "TOX", "CTLA4",
  "CD274", "VSIR", "ENTPD1", "CD96",
  "ARG1", "IL10", "TGFB1",
  "CD163", "MRC1",
  "CXCL8", "CCL2", "VEGFA",
  "SPP1", "MIF"
)

geneset_of_interest <- intersect(geneset_of_interest_all, expressed_genes_receiver)
message(sprintf("Target gene set: %d / 20 present in receiver cells (manuscript: 18)",
                length(geneset_of_interest)))

# Background = all expressed genes in receiver
background_expressed_genes <- expressed_genes_receiver

# -----------------------------------------------------------------------------
# 5. Filter NicheNet ligands to those expressed in sender cells
# -----------------------------------------------------------------------------
ligands          <- lr_network$from %>% unique()
expressed_ligands <- intersect(ligands, expressed_genes_sender)
message(sprintf("Expressed ligands in GlycoHigh cells: %d (manuscript: 325)",
                length(expressed_ligands)))

# Filter receptor side
receptors         <- lr_network$to %>% unique()
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

# Keep LR pairs where both expressed
lr_network_filtered <- lr_network %>%
  filter(from %in% expressed_ligands, to %in% expressed_receptors)

# Restrict ligand-target matrix to expressed ligands
ligand_target_matrix_filtered <- ligand_target_matrix[
  rownames(ligand_target_matrix) %in% background_expressed_genes,
  colnames(ligand_target_matrix) %in% expressed_ligands,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# 6. Predict ligand activity (corrected AUPR as ranking metric)
# -----------------------------------------------------------------------------
message("Computing ligand activities...")

ligand_activities <- predict_ligand_activities(
  geneset              = geneset_of_interest,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix_filtered,
  potential_ligands    = expressed_ligands
)

# Sort by corrected AUPR (aupr_corrected) per manuscript
ligand_activities <- ligand_activities %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = row_number())

# Save full results
dir.create("results", showWarnings = FALSE)
write.csv(ligand_activities,
          "results/S22_NicheNet_ligand_activities.csv", row.names = FALSE)

# Print manuscript focal ligands
mif_row  <- ligand_activities %>% filter(test_ligand == "MIF")
spp1_row <- ligand_activities %>% filter(test_ligand == "SPP1")

if (nrow(mif_row) > 0) {
  message(sprintf("MIF  — Rank %d / %d | AUPR = %.4f (manuscript: rank 37, AUPR = 0.091)",
                  mif_row$rank, nrow(ligand_activities), mif_row$aupr_corrected))
}
if (nrow(spp1_row) > 0) {
  message(sprintf("SPP1 — Rank %d / %d | AUPR = %.4f (manuscript: rank 250, AUPR = 0.026)",
                  spp1_row$rank, nrow(ligand_activities), spp1_row$aupr_corrected))
}

# Top ligands (top 30 by aupr_corrected)
top_ligands <- ligand_activities %>%
  top_n(30, aupr_corrected) %>%
  pull(test_ligand)

# -----------------------------------------------------------------------------
# 7. Compute ligand–target regulatory potential for top ligands
# -----------------------------------------------------------------------------
message("Computing ligand-target links for top ligands...")

active_ligand_target_links_df <- lapply(top_ligands, function(lig) {
  get_weighted_ligand_target_links(
    ligand               = lig,
    geneset              = geneset_of_interest,
    ligand_target_matrix = ligand_target_matrix_filtered,
    n                    = 250   # top 250 targets per ligand
  )
}) %>% bind_rows() %>% drop_na()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df     = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix_filtered,
  cutoff               = 0.33
)

# Ligand–receptor links
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  ligands              = top_ligands,
  expressed_receptors  = expressed_receptors,
  lr_network           = lr_network_filtered,
  weighted_networks    = weighted_networks$lr_sig
)

# -----------------------------------------------------------------------------
# 8. Generate Supplementary Figure S22
# -----------------------------------------------------------------------------
message("Generating Supplementary Figure S22...")

# ── Panel A: Top ligands ranked by corrected AUPR ──
n_show   <- 30
plot_df  <- ligand_activities %>%
  slice_max(aupr_corrected, n = n_show) %>%
  mutate(
    test_ligand = factor(test_ligand, levels = rev(test_ligand)),
    highlight = case_when(
      test_ligand == "MIF"  ~ "MIF (rank 37)",
      test_ligand == "SPP1" ~ "SPP1 (rank 250)",
      TRUE                  ~ "Other"
    )
  )

p_activity <- ggplot(plot_df,
                     aes(x = aupr_corrected, y = test_ligand, fill = highlight)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_vline(xintercept = mean(ligand_activities$aupr_corrected),
             linetype = "dashed", color = "grey40", linewidth = 0.5) +
  scale_fill_manual(
    values = c("MIF (rank 37)" = "#5B4EA6",
               "SPP1 (rank 250)" = "#E07B00",
               "Other" = "#90B8D4"),
    name = "Ligand"
  ) +
  labs(
    title = "Top 30 Predicted Ligands from GlycoHigh Hepatocytes",
    subtitle = paste0("Sender: GlycoHigh hepatocytes (n = ", length(glyco_high_cells),
                      "); Receiver: T/NK + Myeloid (n = ", length(receiver_cells), ")"),
    x = "Corrected AUPR (ligand activity)",
    y = "Ligand"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 11),
    axis.text.y = element_text(
      face = ifelse(levels(plot_df$test_ligand) %in% c("MIF", "SPP1"), "bold", "plain"),
      color = ifelse(levels(plot_df$test_ligand) == "MIF", "#5B4EA6",
               ifelse(levels(plot_df$test_ligand) == "SPP1", "#E07B00", "black"))
    )
  )

# ── Panel B: MIF and SPP1 context relative to all expressed ligands ──
rank_summary <- ligand_activities %>%
  mutate(percentile = 1 - (rank - 1) / max(rank)) %>%
  filter(test_ligand %in% c("MIF", "SPP1"))

p_context <- ggplot(ligand_activities, aes(x = rank, y = aupr_corrected)) +
  geom_line(color = "grey70", linewidth = 0.6) +
  geom_point(data = ligand_activities %>% filter(test_ligand == "MIF"),
             aes(color = "MIF"), size = 3) +
  geom_point(data = ligand_activities %>% filter(test_ligand == "SPP1"),
             aes(color = "SPP1"), size = 3) +
  ggrepel::geom_label_repel(
    data = ligand_activities %>% filter(test_ligand %in% c("MIF", "SPP1")),
    aes(label = sprintf("%s\nRank %d\nAUPR=%.3f", test_ligand, rank, aupr_corrected),
        color = test_ligand),
    size = 3, show.legend = FALSE, box.padding = 0.5
  ) +
  scale_color_manual(values = c("MIF" = "#5B4EA6", "SPP1" = "#E07B00"), name = "Ligand") +
  labs(
    title = "MIF and SPP1 Ranked Among All Expressed Ligands",
    x = "Ligand Rank (by corrected AUPR)",
    y = "Corrected AUPR"
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11),
        legend.position = "bottom")

# ── Panel C: Heatmap of top ligand–target links ──
# Keep only ligands with at least one target link
top_ligands_filtered <- intersect(
  top_ligands,
  unique(active_ligand_target_links_df$ligand)
)

if (length(top_ligands_filtered) > 0 && nrow(active_ligand_target_links) > 0) {
  # Restrict heatmap to top 20 for readability
  ht_ligands <- top_ligands_filtered[seq_len(min(20, length(top_ligands_filtered)))]
  ht_mat <- active_ligand_target_links[
    rownames(active_ligand_target_links) %in% geneset_of_interest,
    colnames(active_ligand_target_links) %in% ht_ligands,
    drop = FALSE
  ]

  ht_df <- as.data.frame(as.table(ht_mat)) %>%
    rename(Target = Var1, Ligand = Var2, Score = Freq) %>%
    filter(!is.na(Score))

  p_heatmap <- ggplot(ht_df, aes(x = Ligand, y = Target, fill = Score)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colors = c("white", "#DAE8F5", "#5B9EC9", "#1A4E7C"),
      name   = "Regulatory\npotential"
    ) +
    labs(
      title = "Ligand–Target Regulatory Potential (Top 20 Ligands)",
      x = "Sender Ligand (GlycoHigh hepatocytes)",
      y = "Target Gene (Receiver immune cells)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x    = element_text(angle = 45, hjust = 1,
                                    face = ifelse(ht_df$Ligand %in% c("MIF","SPP1") %>% unique(),
                                                  "bold", "plain")),
      plot.title     = element_text(face = "bold", size = 11),
      legend.position = "right"
    )
} else {
  p_heatmap <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "No ligand-target links to display") +
    theme_void()
}

# ── Assemble figure ──
fig_s22 <- (p_activity | p_context) / p_heatmap +
  plot_annotation(
    title = paste0(
      "Supplementary Figure S22. NicheNet cross-method validation of immunosuppressive ",
      "ligand activity from GlycoHigh hepatocytes"
    ),
    caption = paste0(
      "NicheNet v2.2.1 (nichenetr). Sender: GlycoHigh hepatocytes (n = ",
      length(glyco_high_cells), "); Receiver: T/NK + Myeloid (n = ",
      length(receiver_cells), ").\n",
      "Target gene set: ", length(geneset_of_interest), " / 20 immunosuppression markers present in receiver cells.\n",
      if (nrow(mif_row)  > 0) sprintf("MIF  rank %d (AUPR = %.3f). ", mif_row$rank,  mif_row$aupr_corrected)  else "",
      if (nrow(spp1_row) > 0) sprintf("SPP1 rank %d (AUPR = %.3f).", spp1_row$rank, spp1_row$aupr_corrected) else "",
      "\nDashed line = mean AUPR across all ", nrow(ligand_activities), " expressed ligands."
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 11),
      plot.caption = element_text(size = 8.5, color = "grey40")
    )
  )

dir.create("figures", showWarnings = FALSE)
ggsave("figures/Supplementary_Figure_S22.pdf", fig_s22,
       width = 16, height = 14, dpi = 300, useDingbats = FALSE)
ggsave("figures/Supplementary_Figure_S22.png", fig_s22,
       width = 16, height = 14, dpi = 300)
message("Supplementary Figure S22 saved to figures/")

# -----------------------------------------------------------------------------
# 9. Save ligand-target and ligand-receptor link tables
# -----------------------------------------------------------------------------
write.csv(active_ligand_target_links_df,
          "results/S22_NicheNet_ligand_target_links.csv", row.names = FALSE)
write.csv(ligand_receptor_links_df,
          "results/S22_NicheNet_ligand_receptor_links.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           "results/15_NicheNet_session_info.txt")
message("\n=== Script 15 complete ===")
