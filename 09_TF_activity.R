# ============================================================
# Script 09: Transcription Factor Activity Analysis
# Section 2.12, Results Section 3.12
# Method: DoRothEA (A+B confidence) + VIPER — 118 TFs
# Generates: Supplementary Figure S16
#
# NOTE: This script replaces the original AUCell-based 5-TF
# analysis. The manuscript describes DoRothEA + VIPER with
# 118 TFs as the final method.
# ============================================================

library(Seurat)
library(dorothea)   # v1.22.0 (Bioconductor)
library(viper)      # v1.44.0 (Bioconductor)
library(ggplot2)
library(dplyr)

setwd("D:/scRNA_project")

# ── Load tumor hepatocytes ────────────────────────────────
seu     <- readRDS("seurat_final.rds")
seu_hep <- subset(seu, cell_type == "Hepatocytes" & site == "Tumor")
cat("Tumor hepatocytes:", ncol(seu_hep), "\n")   # expect 15,391

if (!"glycolysis_score" %in% colnames(seu_hep@meta.data)) {
  stop("glycolysis_score not found. Run 02_glycolysis_scoring.R first.")
}
glyco_scores <- seu_hep$glycolysis_score

# ── Load DoRothEA regulons (A + B confidence) ─────────────
data(dorothea_hs, package = "dorothea")
regulons_ab <- dorothea_hs %>% filter(confidence %in% c("A", "B"))
cat("TF-target interactions (A+B):", nrow(regulons_ab), "\n")

# Convert to VIPER format
viper_regulons <- lapply(split(regulons_ab, regulons_ab$tf), function(x) {
  list(tfmode     = setNames(x$mor, x$target),
       likelihood = rep(1, nrow(x)))
})

# ── Run VIPER ─────────────────────────────────────────────
cat("Running VIPER...\n")
expr_mat <- as.matrix(GetAssayData(seu_hep, assay = "RNA", layer = "data"))

tf_activity <- viper(
  eset        = expr_mat,
  regulon     = viper_regulons,
  eset.filter = FALSE,
  minsize     = 4,
  verbose     = FALSE
)
cat("TFs retained (minsize=4):", nrow(tf_activity), "\n")  # expect ~118

# ── Spearman correlations ─────────────────────────────────
cor_results <- data.frame(TF = rownames(tf_activity),
                          rho = NA_real_, pval = NA_real_)
for (i in seq_len(nrow(tf_activity))) {
  ct <- cor.test(tf_activity[i, ], glyco_scores,
                 method = "spearman", exact = FALSE)
  cor_results$rho[i]  <- ct$estimate
  cor_results$pval[i] <- ct$p.value
}
cor_results$padj <- p.adjust(cor_results$pval, method = "BH")
cor_results      <- cor_results[order(-abs(cor_results$rho)), ]

cat("Significant TFs (padj<0.05):", sum(cor_results$padj < 0.05), "\n")

# Manuscript cross-check
for (tf in c("HIF1A","MYC","FOSL1","FOS","HNF1A","SP1")) {
  r <- cor_results[cor_results$TF == tf, ]
  if (nrow(r) > 0)
    cat(sprintf("  %s rho=%.3f padj=%.2e rank=%d\n",
                tf, r$rho, r$padj, which(cor_results$TF == tf)))
}
# Expected: HIF1A~0.530, MYC~0.427, HNF1A~-0.311, SP1 rank~25

write.csv(cor_results, "FigS16_TF_activity_dorothea.csv", row.names = FALSE)

# ── Supplementary Figure S16: Top 20 TFs ─────────────────
top20 <- cor_results %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = 20) %>%
  mutate(TF        = factor(TF, levels = rev(TF)),
         Direction = ifelse(rho > 0, "Positive", "Negative"),
         sig_label = ifelse(padj < 0.05, "*", ""))

p_s16 <- ggplot(top20, aes(x = rho, y = TF, fill = Direction)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_text(aes(label = sig_label,
                x = rho + sign(rho) * 0.005),
            size = 5, hjust = 0) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  scale_fill_manual(values = c("Positive" = "#E64B35",
                                "Negative" = "#4DBBD5")) +
  labs(
    x        = "Spearman rho (TF activity vs. glycolysis score)",
    y        = NULL,
    title    = "Genome-wide TF Activity Screen (DoRothEA A+B, VIPER)",
    subtitle = paste0("Top 20 of ", nrow(cor_results), " TFs | ",
                      "Tumor hepatocytes n=", ncol(seu_hep),
                      " | * padj<0.05")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom",
        plot.title      = element_text(face = "bold"))

ggsave("FigS16_TF_activity_DoRothEA.pdf", p_s16, width = 8, height = 7)
ggsave("FigS16_TF_activity_DoRothEA.png", p_s16, width = 8, height = 7, dpi = 300)
cat("Supplementary Figure S16 saved.\n\n=== Script 09 complete ===\n")
