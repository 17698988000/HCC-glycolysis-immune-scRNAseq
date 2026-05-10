# =============================================================================
# 14_OXPHOS_metabolic_specificity.R
# Metabolic Pathway Specificity Analysis (Section 2.15)
# Computes OXPHOS AUCell scores and partial Spearman correlations between
# glycolysis/OXPHOS and immunosuppressive ligands (SPP1, MIF) in tumor hepatocytes
# Generates: Supplementary Figure S22
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(ppcor)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Load processed data (tumor-derived hepatocytes from 01_QC_clustering.R)
# -----------------------------------------------------------------------------
message("Loading tumor hepatocyte data...")

# Assumes the full Seurat object has been saved by 01_QC_clustering.R
seurat_obj <- readRDS("results/seurat_hcc_annotated.rds")

# Subset to tumor-derived hepatocytes (n = 15,391 as in manuscript)
tumor_hep <- subset(seurat_obj,
                    subset = cell_type == "Hepatocyte" & tissue == "Tumor")
message(sprintf("Tumor hepatocytes retained: %d", ncol(tumor_hep)))

# -----------------------------------------------------------------------------
# 2. Define OXPHOS gene set (39 core subunit-encoding genes, Section 2.15)
# -----------------------------------------------------------------------------
oxphos_genes <- c(
  # Complex I (NADH dehydrogenase)
  "NDUFB3", "NDUFB4", "NDUFB5", "NDUFB6", "NDUFB7", "NDUFB8", "NDUFB9",
  "NDUFB10", "NDUFA1", "NDUFA2", "NDUFA3", "NDUFA4",
  "NDUFS1", "NDUFS2", "NDUFS3",
  # Complex II (succinate dehydrogenase)
  "SDHA", "SDHB", "SDHC", "SDHD",
  # Complex III (cytochrome bc1)
  "UQCRB", "UQCRC1", "UQCRC2", "UQCRQ", "CYC1",
  # Complex IV (cytochrome c oxidase)
  "COX4I1", "COX5A", "COX5B", "COX6A1", "COX6B1", "COX7A2", "COX8A",
  # Complex V (ATP synthase)
  "ATP5F1A", "ATP5F1B", "ATP5F1C", "ATP5F1D", "ATP5F1E",
  "ATP5MC1", "ATP5MC2", "ATP5MC3", "ATP5PB"
)

# Glycolysis gene set (22 genes, same as Section 2.4 / 02_glycolysis_scoring.R)
glycolysis_genes <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKM",
  "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1",
  "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
  "SLC2A1", "SLC2A3", "PFKFB3", "GCK"
)

# Filter to genes present in the dataset
oxphos_genes     <- intersect(oxphos_genes,     rownames(tumor_hep))
glycolysis_genes <- intersect(glycolysis_genes, rownames(tumor_hep))
message(sprintf("OXPHOS genes used: %d / 39", length(oxphos_genes)))
message(sprintf("Glycolysis genes used: %d / 22", length(glycolysis_genes)))

# -----------------------------------------------------------------------------
# 3. Compute AUCell scores for OXPHOS
# (Glycolysis AUCell scores assumed already in tumor_hep metadata from Script 02)
# -----------------------------------------------------------------------------
message("Computing OXPHOS AUCell scores...")
expr_matrix <- GetAssayData(tumor_hep, slot = "counts")

# Build gene sets
gene_sets <- list(
  Glycolysis = glycolysis_genes,
  OXPHOS     = oxphos_genes
)

# Run AUCell
cells_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)
cells_AUC      <- AUCell_calcAUC(gene_sets, cells_rankings,
                                  aucMaxRank = ceiling(0.05 * nrow(expr_matrix)),
                                  verbose = FALSE)

# Extract scores
auc_df <- as.data.frame(t(getAUC(cells_AUC)))
colnames(auc_df) <- c("Glycolysis_AUC", "OXPHOS_AUC")

# Merge with existing metadata (use glycolysis score from Script 02 if available,
# otherwise use freshly computed one)
if ("Glycolysis_AUC" %in% colnames(tumor_hep@meta.data)) {
  auc_df$Glycolysis_AUC <- tumor_hep@meta.data[rownames(auc_df), "Glycolysis_AUC"]
  message("Using pre-computed glycolysis AUC scores from metadata.")
}

tumor_hep <- AddMetaData(tumor_hep, metadata = auc_df)

# Extract expression of immunosuppressive ligands
ligands <- c("SPP1", "MIF")
ligand_expr <- as.data.frame(
  t(GetAssayData(tumor_hep, slot = "data")[ligands, , drop = FALSE])
)

# Combine into analysis dataframe
analysis_df <- cbind(auc_df, ligand_expr)
message(sprintf("Analysis dataframe: %d cells x %d variables", nrow(analysis_df), ncol(analysis_df)))

# -----------------------------------------------------------------------------
# 4. Check independence of glycolysis and OXPHOS programs
# -----------------------------------------------------------------------------
message("Checking glycolysis-OXPHOS correlation...")
cor_go <- cor.test(analysis_df$Glycolysis_AUC, analysis_df$OXPHOS_AUC,
                   method = "spearman", exact = FALSE)
message(sprintf("Glycolysis vs OXPHOS: rho = %.3f, p = %.2e",
                cor_go$estimate, cor_go$p.value))
# Should be ~0.066 per manuscript

# -----------------------------------------------------------------------------
# 5. Partial Spearman correlations (ppcor package)
# For each ligand:
#   (a) corr(ligand ~ Glycolysis | OXPHOS)
#   (b) corr(ligand ~ OXPHOS | Glycolysis)
# -----------------------------------------------------------------------------
message("Computing partial Spearman correlations...")

results <- list()

for (ligand in ligands) {
  dat <- data.frame(
    ligand    = analysis_df[[ligand]],
    Glyco     = analysis_df$Glycolysis_AUC,
    OXPHOS    = analysis_df$OXPHOS_AUC
  )

  # Partial correlation: ligand ~ Glycolysis controlling for OXPHOS
  pc_glyco <- pcor.test(dat$ligand, dat$Glyco, dat$OXPHOS, method = "spearman")

  # Partial correlation: ligand ~ OXPHOS controlling for Glycolysis
  pc_oxphos <- pcor.test(dat$ligand, dat$OXPHOS, dat$Glyco, method = "spearman")

  results[[ligand]] <- data.frame(
    Ligand        = ligand,
    Variable      = c("Glycolysis (ctrl OXPHOS)", "OXPHOS (ctrl Glycolysis)"),
    Partial_rho   = c(pc_glyco$estimate,  pc_oxphos$estimate),
    p_value       = c(pc_glyco$p.value,   pc_oxphos$p.value),
    Statistic     = c(pc_glyco$statistic, pc_oxphos$statistic),
    n             = c(pc_glyco$n,         pc_oxphos$n)
  )

  message(sprintf(
    "%s — Glycolysis partial rho = %.3f (p = %.2e) | OXPHOS partial rho = %.3f (p = %.2e)",
    ligand,
    pc_glyco$estimate,  pc_glyco$p.value,
    pc_oxphos$estimate, pc_oxphos$p.value
  ))
}

results_df <- do.call(rbind, results)
results_df$padj <- p.adjust(results_df$p_value, method = "BH")
results_df$neg_log10_p <- -log10(results_df$p_value)

# Save table
dir.create("results", showWarnings = FALSE)
write.csv(results_df, "results/S22_partial_correlations.csv", row.names = FALSE)
message("Partial correlation results saved to results/S22_partial_correlations.csv")

# Print summary matching manuscript values
message("\n=== Manuscript cross-check ===")
message("SPP1 ~ Glycolysis | OXPHOS  expected rho ≈  0.227, p ≈ 3.43e-131")
message("SPP1 ~ OXPHOS | Glycolysis  expected rho ≈  0.125, p ≈ 6.44e-40")
message("MIF  ~ Glycolysis | OXPHOS  expected rho ≈  0.105, p ≈ 6.72e-29")
message("MIF  ~ OXPHOS | Glycolysis  expected rho ≈  0.242, p ≈ 2.70e-148")

# -----------------------------------------------------------------------------
# 6. Generate Supplementary Figure S22
# -----------------------------------------------------------------------------
message("Generating Supplementary Figure S22...")

# Panel A: Glycolysis vs OXPHOS scatter (confirm independence)
p_independence <- ggplot(analysis_df[sample(nrow(analysis_df), min(3000, nrow(analysis_df))), ],
                         aes(x = Glycolysis_AUC, y = OXPHOS_AUC)) +
  geom_point(alpha = 0.3, size = 0.5, color = "#555555") +
  geom_smooth(method = "lm", color = "#E04040", linewidth = 0.8) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = sprintf("Spearman rho = %.3f\np = %.2e",
                           cor_go$estimate, cor_go$p.value),
           size = 3.5) +
  labs(title = "Glycolysis vs. OXPHOS in Tumor Hepatocytes",
       subtitle = sprintf("n = %d cells", nrow(analysis_df)),
       x = "Glycolysis AUCell Score",
       y = "OXPHOS AUCell Score") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 11))

# Panel B: Partial correlation bar plot (main figure)
plot_df <- results_df %>%
  mutate(
    Direction = ifelse(Partial_rho > 0, "Positive", "Negative"),
    Label = sprintf("rho = %.3f\np = %.2e", Partial_rho, p_value),
    Metabolic_pathway = ifelse(grepl("Glycolysis", Variable), "Glycolysis\n(ctrl OXPHOS)", "OXPHOS\n(ctrl Glycolysis)")
  )

# Color scheme: SPP1 = orange/teal, MIF = purple/green
p_partial <- ggplot(plot_df, aes(x = Metabolic_pathway, y = Partial_rho,
                                  fill = Ligand)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "white") +
  geom_text(aes(label = sprintf("rho=%.3f", Partial_rho),
                y = Partial_rho + sign(Partial_rho) * 0.005),
            position = position_dodge(width = 0.8),
            size = 3, vjust = -0.3) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c("SPP1" = "#E07B00", "MIF" = "#5B4EA6"),
                    name = "Ligand") +
  facet_wrap(~Ligand, ncol = 2) +
  labs(
    title = "Partial Spearman Correlations: Metabolic Programs vs. Immunosuppressive Ligands",
    subtitle = "Tumor-derived hepatocytes (n = 15,391)",
    x = "Metabolic Pathway (covariate controlled)",
    y = "Partial Spearman rho"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    strip.text    = element_text(face = "bold"),
    legend.position = "none"
  )

# Panel C: Signed -log10(p) to emphasize significance
p_significance <- ggplot(plot_df,
                         aes(x = Metabolic_pathway,
                             y = sign(Partial_rho) * neg_log10_p,
                             fill = Ligand)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "white") +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  geom_hline(yintercept = c(-log10(0.05), log10(0.05)),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  scale_fill_manual(values = c("SPP1" = "#E07B00", "MIF" = "#5B4EA6")) +
  facet_wrap(~Ligand, ncol = 2) +
  labs(
    title = "Statistical Significance of Partial Correlations",
    x = "Metabolic Pathway",
    y = "sign(rho) × -log10(p)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold", size = 11),
    strip.text    = element_text(face = "bold"),
    legend.position = "none"
  )

# Combine panels
fig_s22 <- (p_independence | p_partial) / p_significance +
  plot_annotation(
    title = "Supplementary Figure S22. Glycolysis and OXPHOS independently associate with distinct immunosuppressive ligands",
    caption = paste0(
      "Partial Spearman correlations computed using ppcor v1.1.\n",
      "SPP1 is predominantly glycolysis-associated; MIF is predominantly OXPHOS-associated.\n",
      "Glycolysis-OXPHOS correlation: rho = ", round(cor_go$estimate, 3),
      ", p = ", formatC(cor_go$p.value, format = "e", digits = 2)
    ),
    theme = theme(
      plot.title   = element_text(face = "bold", size = 12),
      plot.caption = element_text(size = 9, color = "grey40")
    )
  )

# Save figure
dir.create("figures", showWarnings = FALSE)
ggsave("figures/Supplementary_Figure_S22.pdf", fig_s22,
       width = 14, height = 10, dpi = 300, useDingbats = FALSE)
ggsave("figures/Supplementary_Figure_S22.png", fig_s22,
       width = 14, height = 10, dpi = 300)
message("Supplementary Figure S22 saved to figures/")

# -----------------------------------------------------------------------------
# 7. Save session info
# -----------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()),
           "results/14_OXPHOS_session_info.txt")
message("\n=== Script 14 complete ===")
