# ============================================================
# 12_drug_repurposing.R
# Purpose: Drug repurposing analysis using Fisher's exact test
#          enrichment against curated drug-target database
# Corresponds to: Methods Section 2.15, Results Section 3.14
# ============================================================

library(limma)
library(ggplot2)

# ── Load TCGA-LIHC expression and risk scores ──────────────
expr <- read.table("D:/scRNA_project/TCGA_LIHC_expression.txt",
                   header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
risk_data <- read.csv("D:/scRNA_project/TCGA_LIHC_risk_scores.csv")

# ── DEG analysis: high vs low risk group ──────────────────
common_samples <- intersect(colnames(expr), risk_data$sample)
expr_sub  <- expr[, common_samples]
risk_sub  <- risk_data[match(common_samples, risk_data$sample), ]

group <- factor(risk_sub$risk_group, levels = c("Low","High"))
design <- model.matrix(~ group)
fit    <- lmFit(expr_sub, design)
fit    <- eBayes(fit)
deg    <- topTable(fit, coef = 2, number = Inf, adjust.method = "BH")
deg$gene <- rownames(deg)

# Save DEG table
write.csv(deg, "D:/scRNA_project/TCGA_LIHC_risk_DEG_corrected.csv")
cat("Total DEGs (FDR<0.05):", sum(deg$adj.P.Val < 0.05), "\n")

# Upregulated in high-risk group
up_genes <- rownames(deg[deg$logFC > 0.5 & deg$adj.P.Val < 0.05, ])
cat("Upregulated genes (logFC>0.5, FDR<0.05):", length(up_genes), "\n")

# ── Curated drug-target database ──────────────────────────
# Sources: DrugBank, Therapeutic Target Database (TTD),
#          published HCC pharmacology literature
drug_db <- list(
  `2-DG`          = c("HK1","HK2","SLC2A1","SLC2A3","PFKP","PFKL","PFKM",
                       "ENO1","LDHA","PKM","ALDOA","GPI","TPI1"),
  Lonidamine      = c("HK2","HK1","LDHA","MCT4","SLC16A3","VDAC1","PKM","ENO1"),
  `SB-265610`     = c("CXCR2","CXCL8","CXCL1","CXCL5","IL8"),
  Marimastat      = c("MMP1","MMP2","MMP3","MMP7","MMP9","MMP14","ADAM10","ADAM17"),
  Reparixin       = c("CXCR1","CXCR2","CXCL8","CXCL1","CXCL2","IL6","TNF"),
  Celecoxib       = c("PTGS2","PTGER4","PTGER2","CXCL8","IL6","SPP1","MIF",
                      "VEGFA","BCL2"),
  Ibuprofen       = c("PTGS1","PTGS2","PTGER4","PTGER2","CXCL8","IL6","TNF",
                      "MIF","SPP1"),
  Metformin       = c("PRKAA1","PRKAA2","MTOR","HIF1A","LDHA","SLC2A1","PKM",
                      "VEGFA","IGF1R"),
  Rapamycin       = c("MTOR","RPTOR","RPS6KB1","EIF4EBP1","HIF1A","LDHA",
                      "SLC2A1","ENO1","PKM"),
  Cabozantinib    = c("MET","VEGFR2","AXL","RET","KIT","FLT3","TRKB","TYRO3",
                      "HGF","SPP1"),
  Sorafenib       = c("VEGFA","VEGFB","PDGFRA","PDGFRB","RAF1","BRAF","KIT",
                      "FLT3","RET","MAPK1","MAPK3","HIF1A","SLC2A1","LDHA",
                      "PKM","ENO1"),
  Bevacizumab     = c("VEGFA","VEGFB","VEGFC","VEGFD","HIF1A","SLC2A1"),
  Tasquinimod     = c("S100A9","ANXA2","HIF1A","VEGFA","HDAC4","MMP9"),
  Pexidartinib    = c("CSF1R","KIT","FLT3","SPP1","CD68","MRC1","ARG1"),
  Cabiralizumab   = c("CSF1R","CD68","MRC1","ARG1","SPP1","IL10","TGFB1"),
  Gefitinib       = c("EGFR","ERBB2","MAPK1","AKT1","MTOR","SLC2A1","HIF1A"),
  Vorinostat      = c("HDAC1","HDAC2","HDAC3","HDAC6","HIF1A","VEGFA",
                      "SLC2A1","LDHA"),
  Galunisertib    = c("TGFBR1","TGFB1","SMAD2","SMAD3","FOXP3","CD8A","IFNG"),
  Tiragolumab     = c("TIGIT","CD226","PVRL2","CD274","CD8A","IFNG"),
  Bintrafusp      = c("TGFB1","TGFB2","TGFB3","PDCD1LG2","CD274","SMAD2","SMAD3"),
  Entinostat      = c("HDAC1","HDAC2","HDAC3","FOXP3","CD8A","IFNG","IL10",
                      "TGFB1","ARG1"),
  Nivolumab       = c("PDCD1","CD274","CD8A","CD8B","HAVCR2","LAG3","TIGIT","IFNG"),
  Pembrolizumab   = c("PDCD1","CD274","CD8A","GZMB","PRF1","IFNG","LAG3"),
  Atezolizumab    = c("CD274","PDCD1LG2","CD8A","CD8B","GZMB","PRF1","IFNG",
                      "TIGIT","LAG3"),
  Tremelimumab    = c("CTLA4","CD80","CD86","FOXP3","IL2","IFNG"),
  Ipilimumab      = c("CTLA4","CD80","CD86","FOXP3","IL2","CD4","ICOS"),
  Lenvatinib      = c("VEGFR1","VEGFR2","VEGFR3","FGFR1","FGFR2","FGFR3",
                      "FGFR4","PDGFRA","KIT","RET","VEGFA","FGF19","HGF"),
  Regorafenib     = c("VEGFR1","VEGFR2","VEGFR3","PDGFRA","PDGFRB","RAF1",
                      "BRAF","KIT","RET","FGFR1","TIE2","MAPK1"),
  Ramucirumab     = c("VEGFR2","VEGFA","VEGFB","VEGFC"),
  Galunisertib2   = c("TGFBR1","TGFB1","SMAD2","SMAD3","VEGFA","MMP9","LDHA")
)

# ── Fisher enrichment test ────────────────────────────────
total_genes <- nrow(deg)

results <- lapply(names(drug_db), function(drug_name) {
  drug_genes <- drug_db[[drug_name]]
  a <- sum(up_genes %in% drug_genes)
  if (a == 0) return(NULL)
  b  <- length(up_genes) - a
  c_ <- sum(drug_genes %in% rownames(deg)) - a
  d  <- total_genes - a - b - c_
  p  <- fisher.test(matrix(c(a, b, c_, d), 2, 2),
                    alternative = "greater")$p.value
  overlap_genes <- paste(up_genes[up_genes %in% drug_genes], collapse = ",")
  data.frame(drug = drug_name, overlap = a, p_value = p,
             overlap_genes = overlap_genes, stringsAsFactors = FALSE)
})

res_df      <- do.call(rbind, Filter(Negate(is.null), results))
res_df$FDR  <- p.adjust(res_df$p_value, method = "BH")
res_df      <- res_df[order(res_df$FDR), ]

cat("Significant drugs (FDR<0.05):", sum(res_df$FDR < 0.05), "\n")
print(res_df[res_df$FDR < 0.05, c("drug","overlap","FDR","overlap_genes")])
write.csv(res_df, "D:/scRNA_project/drug_repurposing_results.csv", row.names = FALSE)

# ── Plot FigS19 ───────────────────────────────────────────
plot_df <- res_df[res_df$FDR < 0.05, ]
plot_df$neg_log_FDR <- -log10(plot_df$FDR)
plot_df <- plot_df[order(plot_df$neg_log_FDR, decreasing = TRUE), ]

# Evidence classification
fda_hcc <- c("Sorafenib","Bevacizumab","Atezolizumab","Nivolumab",
             "Pembrolizumab","Lenvatinib","Regorafenib","Ramucirumab","Cabozantinib")
plot_df$category <- "Clinical/Investigational"
plot_df$category[plot_df$drug %in% fda_hcc] <- "FDA-approved (HCC)"
plot_df$category[plot_df$drug %in% c("2-DG","Lonidamine")] <- "Preclinical"

plot_df$drug <- factor(plot_df$drug, levels = plot_df$drug)

color_map <- c("FDA-approved (HCC)"  = "firebrick3",
               "Clinical/Investigational" = "steelblue3",
               "Preclinical"         = "grey50")

p <- ggplot(plot_df, aes(x = neg_log_FDR, y = drug, fill = category)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", linewidth = 0.3) +
  geom_text(aes(label = overlap_genes),
            x = 0.05, hjust = 0, size = 2.8, color = "white", fontface = "bold") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "black", linewidth = 0.5) +
  scale_fill_manual(values = color_map, name = "Evidence level") +
  scale_x_continuous(expand = c(0, 0),
                     limits = c(0, max(plot_df$neg_log_FDR) * 1.15)) +
  labs(x = expression(-log[10](FDR)), y = NULL,
       title = "Drug Repurposing: Glycolysis-High HCC Signature",
       subtitle = paste0("Fisher's exact test | n = ", length(up_genes),
                         " upregulated genes | ",
                         sum(res_df$FDR < 0.05), " significant compounds (FDR<0.05)")) +
  theme_classic(base_size = 12) +
  theme(axis.text.y   = element_text(size = 11, color = "black"),
        legend.position = c(0.75, 0.25),
        plot.title    = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9, color = "grey40"))

ggsave("D:/scRNA_project/FigS19_drug_repurposing_v3.png",
       p, width = 9, height = 7, dpi = 300)
ggsave("D:/scRNA_project/FigS19_drug_repurposing_v3.pdf",
       p, width = 9, height = 7)

cat("Drug repurposing analysis complete.\n")
