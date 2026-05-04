# ============================================================
# 11_spatial_transcriptomics.R
# Purpose: Spatial validation of glycolysis-immunosuppression
#          co-localization in HCC Visium samples (GSE238264)
# Corresponds to: Methods Section 2.14, Results Section 3.13
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# ── Glycolysis gene set ────────────────────────────────────
glyco_genes <- c("HK1","HK2","GPI","PFKL","PFKP","PFKM","ALDOA","ALDOB","ALDOC",
                 "TPI1","GAPDH","PGK1","PGAM1","ENO1","ENO2","PKM","LDHA","LDHB",
                 "SLC2A1","SLC2A3","PFKFB3","GCK")

# ── Sample IDs ────────────────────────────────────────────
sample_ids <- c("HCC1R","HCC2R","HCC3R","HCC4R")

# NOTE: Download GSE238264 data from GEO and set path below
# Each sample folder should contain: barcodes.tsv.gz, features.tsv.gz,
# matrix.mtx.gz, tissue_positions.csv, scalefactors_json.json,
# tissue_hires_image.png
data_dir <- "D:/scRNA_project/GSE238264/"

results_list <- list()

for (sid in sample_ids) {
  cat("Processing sample:", sid, "\n")
  sample_path <- file.path(data_dir, sid)

  # Load Visium data
  seu_st <- Load10X_Spatial(data.dir = sample_path, slice = sid)

  # QC filter
  seu_st <- subset(seu_st,
                   nFeature_Spatial >= 200 &
                   PercentageFeatureSet(seu_st, pattern = "^MT-") < 25)
  cat("  Spots after QC:", ncol(seu_st), "\n")

  # Normalize
  seu_st <- SCTransform(seu_st, assay = "Spatial", verbose = FALSE)

  # Glycolysis module score
  avail_genes <- glyco_genes[glyco_genes %in% rownames(seu_st)]
  seu_st <- AddModuleScore(
    seu_st,
    features = list(avail_genes),
    name     = "Glycolysis_Score",
    ctrl     = 100
  )

  # Stratify spots
  med_score <- median(seu_st$Glycolysis_Score1)
  seu_st$glyco_group <- ifelse(seu_st$Glycolysis_Score1 >= med_score,
                                "GlycoHigh", "GlycoLow")

  # Differential expression of ENO1, SPP1, MIF
  target_genes <- c("ENO1","SPP1","MIF")
  avail_targets <- target_genes[target_genes %in% rownames(seu_st)]

  expr_mat  <- GetAssayData(seu_st, assay = "SCT", layer = "data")
  high_cells <- colnames(seu_st)[seu_st$glyco_group == "GlycoHigh"]
  low_cells  <- colnames(seu_st)[seu_st$glyco_group == "GlycoLow"]

  de_results <- sapply(avail_targets, function(g) {
    wt <- wilcox.test(
      expr_mat[g, high_cells],
      expr_mat[g, low_cells],
      alternative = "greater"
    )
    fc <- mean(expr_mat[g, high_cells]) - mean(expr_mat[g, low_cells])
    c(p_value = wt$p.value, log2FC = fc)
  })

  results_list[[sid]] <- as.data.frame(t(de_results))
  results_list[[sid]]$gene   <- rownames(results_list[[sid]])
  results_list[[sid]]$sample <- sid

  cat("  DE results:\n"); print(results_list[[sid]])

  # Spatial feature plots (Figure 8A)
  p1 <- SpatialFeaturePlot(seu_st, features = "ENO1") + ggtitle(paste(sid, "ENO1"))
  p2 <- SpatialFeaturePlot(seu_st, features = "Glycolysis_Score1") +
    ggtitle(paste(sid, "Glycolysis Score"))
  p3 <- SpatialFeaturePlot(seu_st, features = "SPP1") + ggtitle(paste(sid, "SPP1"))
  p4 <- SpatialFeaturePlot(seu_st, features = "MIF")  + ggtitle(paste(sid, "MIF"))

  p_combined <- (p1 | p2) / (p3 | p4)
  ggsave(paste0("D:/scRNA_project/Fig8A_spatial_", sid, ".png"),
         p_combined, width = 10, height = 8, dpi = 200)
}

# ── Summary bubble plot (Figure 8B / FigS18) ──────────────
all_results <- do.call(rbind, results_list)
all_results$neg_log_p <- -log10(all_results$p_value)
all_results$sig <- ifelse(all_results$p_value < 0.05,
                           ifelse(all_results$p_value < 0.001, "***",
                                  ifelse(all_results$p_value < 0.01, "**", "*")),
                           "ns")

p_bubble <- ggplot(all_results,
                   aes(x = sample, y = gene,
                       size = neg_log_p, color = log2FC)) +
  geom_point() +
  geom_text(aes(label = sig), vjust = -1.2, size = 3.5) +
  scale_size_continuous(name = "-log10(p)", range = c(2, 12)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "log2FC") +
  labs(x = NULL, y = NULL,
       title = "Spatial Validation: Glycolysis-High vs Low Spots",
       subtitle = "GlycoHigh vs GlycoLow | 4 HCC Visium samples") +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(size = 12))

ggsave("D:/scRNA_project/Fig8B_spatial_summary.png",
       p_bubble, width = 8, height = 5, dpi = 300)
ggsave("D:/scRNA_project/Fig8B_spatial_summary.pdf",
       p_bubble, width = 8, height = 5)

# Save results table
write.csv(all_results, "D:/scRNA_project/spatial_DE_results.csv", row.names = FALSE)
cat("Spatial transcriptomics analysis complete.\n")
