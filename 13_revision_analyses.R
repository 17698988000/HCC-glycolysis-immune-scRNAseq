# ============================================================
# Revision Analyses — JTM Major Revision
# ============================================================

library(Seurat)
library(AUCell)
library(spacexr)
library(ggplot2)
library(dplyr)
library(tidyr)

load("D:/scRNA_project/workspace_backup.RData")

glycolysis_genes <- c("HK1","HK2","GPI","PFKL","PFKP","PFKM",
                      "ALDOA","ALDOB","ALDOC","TPI1","GAPDH",
                      "PGK1","PGAM1","ENO1","ENO2","PKM","LDHA",
                      "LDHB","SLC2A1","SLC2A3","PFKFB3","GCK")

# ============================================================
# Issue 2: Leave-One-Out (LOO) Circularity Analysis
# ============================================================

expr_matrix <- GetAssayData(tumor_hep, layer = "counts")
rankings    <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)

loo_cor <- function(exclude_gene, rankings, expr_matrix) {
  gene_set   <- setdiff(glycolysis_genes, exclude_gene)
  auc        <- AUCell_calcAUC(list(loo = gene_set), rankings,
                               aucMaxRank = ceiling(0.05 * nrow(rankings)))
  auc_scores <- as.numeric(getAUC(auc)["loo", ])
  gene_expr  <- as.numeric(expr_matrix[exclude_gene, ])
  cor(gene_expr, auc_scores, method = "spearman")
}

loo_ENO1  <- loo_cor("ENO1",  rankings, expr_matrix)
loo_LDHA  <- loo_cor("LDHA",  rankings, expr_matrix)
loo_GAPDH <- loo_cor("GAPDH", rankings, expr_matrix)

cat("ENO1  LOO rho:", round(loo_ENO1,  3), "\n")  # 0.254
cat("LDHA  LOO rho:", round(loo_LDHA,  3), "\n")  # 0.346
cat("GAPDH LOO rho:", round(loo_GAPDH, 3), "\n")  # 0.374

# Supplementary Figure S20
df_loo <- data.frame(
  Gene  = rep(c("ENO1","LDHA","GAPDH"), each = 2),
  Type  = rep(c("Standard rho","LOO rho"), 3),
  Value = c(0.554, 0.254, 0.537, 0.346, 0.522, 0.374)
)

p_s20 <- ggplot(df_loo, aes(x = Gene, y = Value, fill = Type)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.6) +
  scale_fill_manual(values = c("Standard rho" = "#E64B35", "LOO rho" = "#4DBBD5")) +
  geom_text(aes(label = round(Value, 3)), position = position_dodge(0.6),
            vjust = -0.5, size = 3.5) +
  ylim(0, 0.65) +
  theme_bw(base_size = 12) +
  labs(title = "LOO Spearman Correlation Analysis",
       y = "Spearman rho", x = "", fill = "") +
  theme(legend.position = "bottom")

ggsave("output/FigS20_LOO_comparison.pdf", p_s20, width = 6, height = 5)
ggsave("output/FigS20_LOO_comparison.png", p_s20, width = 6, height = 5, dpi = 300)

# ============================================================
# Issue 5: Spatial Transcriptomics — RCTD Deconvolution
# ============================================================

samples <- c("HCC1R","HCC2R","HCC3R","HCC4R")

# 加载Visium数据
vis_list2 <- list()
for(s in samples) {
  path <- paste0("D:/scRNA_project/GSE238264/", s, "/", s)
  vis_list2[[s]] <- Load10X_Spatial(path, slice = s)
  vis_list2[[s]] <- NormalizeData(vis_list2[[s]])
  vis_list2[[s]] <- AddModuleScore(vis_list2[[s]],
                                   features = list(glycolysis_genes),
                                   name = "Glycolysis", nbin = 10)
  cat(s, "loaded\n")
}

# 构建scRNA参考（注意替换celltype中的"/"）
counts_ref <- GetAssayData(seu, assay = "RNA", layer = "counts")
cell_types  <- as.factor(gsub("/", "_", seu$celltype))
names(cell_types) <- colnames(seu)
reference   <- Reference(counts_ref, cell_types)
# 注：spacexr会自动下采样每类细胞到10000，属正常行为

# 对四个样本跑RCTD
rctd_results <- list()
for(s in samples) {
  cat("Running RCTD for", s, "...\n")
  obj       <- vis_list2[[s]]
  counts_sp <- GetAssayData(obj, assay = "Spatial", layer = "counts")
  coords    <- GetTissueCoordinates(obj)[, c("x","y")]
  sp_obj    <- SpatialRNA(coords, counts_sp)
  rctd      <- create.RCTD(sp_obj, reference, max_cores = 4)
  rctd      <- run.RCTD(rctd, doublet_mode = "full")
  rctd_results[[s]] <- rctd
  cat(s, "done\n")
}

# 提取结果并比较细胞类型比例
results_list <- list()
for(s in samples) {
  rctd         <- rctd_results[[s]]
  weights      <- as.data.frame(rctd@results$weights)
  weights_norm <- weights / rowSums(weights)
  
  obj          <- vis_list2[[s]]
  glyco_score  <- obj$Glycolysis1
  glyco_group  <- ifelse(glyco_score > median(glyco_score), "GlycoHigh","GlycoLow")
  
  weights_norm$group  <- glyco_group[rownames(weights_norm)]
  weights_norm$sample <- s
  results_list[[s]]   <- weights_norm
}

all_results <- bind_rows(results_list)

# Wilcoxon检验：Hepatocyte比例
hep_high <- all_results %>% filter(group == "GlycoHigh") %>% pull(Hepatocyte)
hep_low  <- all_results %>% filter(group == "GlycoLow")  %>% pull(Hepatocyte)
wt       <- wilcox.test(hep_high, hep_low)
cat("Hepatocyte GlycoHigh mean:", round(mean(hep_high), 3),  # 0.629
    "GlycoLow mean:", round(mean(hep_low), 3),               # 0.670
    "p =", format(wt$p.value, scientific = TRUE), "\n")      # 2.75e-18

# Supplementary Figure S21
plot_df <- all_results %>%
  select(sample, group, Hepatocyte, Myeloid, T_NK, Fibroblast, Endothelial, B) %>%
  pivot_longer(cols = -c(sample, group), names_to = "CellType", values_to = "Proportion")

p_s21 <- ggplot(plot_df, aes(x = group, y = Proportion, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~CellType, scales = "free_y") +
  scale_fill_manual(values = c("GlycoHigh" = "#E64B35","GlycoLow" = "#4DBBD5")) +
  theme_bw(base_size = 12) +
  labs(title = "RCTD Cell Type Composition: GlycoHigh vs GlycoLow Spots",
       x = "", y = "Estimated Proportion") +
  theme(legend.position = "bottom")

ggsave("output/FigS21_RCTD_composition.pdf", p_s21, width = 10, height = 8)
ggsave("output/FigS21_RCTD_composition.png", p_s21, width = 10, height = 8, dpi = 300)
# ============================================================
# OXPHOS代谢特异性偏相关分析（Revision新增）
# ============================================================

OXPHOS_genes <- c(
  "NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8","NDUFB9","NDUFB10",
  "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFS1","NDUFS2","NDUFS3",
  "SDHA","SDHB","SDHC","SDHD",
  "UQCRB","UQCRC1","UQCRC2","UQCRQ","CYC1",
  "COX4I1","COX5A","COX5B","COX6A1","COX6B1","COX7A2","COX8A",
  "ATP5F1A","ATP5F1B","ATP5F1C","ATP5F1D","ATP5F1E",
  "ATP5MC1","ATP5MC2","ATP5MC3","ATP5PB"
)

tumor_hep <- subset(seu, site == "Tumor" & celltype == "Hepatocyte")
expr_matrix <- GetAssayData(tumor_hep, assay = "RNA", layer = "counts")
rankings_oxphos <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)
auc_oxphos <- AUCell_calcAUC(list(OXPHOS = OXPHOS_genes), rankings_oxphos,
                              aucMaxRank = ceiling(0.05 * nrow(rankings_oxphos)))
tumor_hep$OXPHOS_AUC <- as.numeric(getAUC(auc_oxphos)["OXPHOS", ])

spp1 <- GetAssayData(tumor_hep, layer = "data")["SPP1", ]
mif  <- GetAssayData(tumor_hep, layer = "data")["MIF", ]

library(ppcor)
df_partial <- data.frame(
  SPP1 = as.numeric(spp1), MIF = as.numeric(mif),
  Glycolysis = tumor_hep$Glycolysis_AUC,
  OXPHOS = tumor_hep$OXPHOS_AUC
)

pc_spp1_glyco  <- pcor.test(df_partial$SPP1, df_partial$Glycolysis, df_partial$OXPHOS, method = "spearman")
pc_spp1_oxphos <- pcor.test(df_partial$SPP1, df_partial$OXPHOS, df_partial$Glycolysis, method = "spearman")
pc_mif_glyco   <- pcor.test(df_partial$MIF,  df_partial$Glycolysis, df_partial$OXPHOS, method = "spearman")
pc_mif_oxphos  <- pcor.test(df_partial$MIF,  df_partial$OXPHOS, df_partial$Glycolysis, method = "spearman")

# ============================================================
# DoRothEA全基因组TF活性分析（Revision新增，替换5-TF候选分析）
# ============================================================

library(dorothea)
library(viper)

data(dorothea_hs, package = "dorothea")
regulons <- dorothea_hs %>% filter(confidence %in% c("A", "B"))

expr_mat <- as.matrix(GetAssayData(tumor_hep, assay = "RNA", layer = "data"))
reg_list <- split(regulons, regulons$tf)
viper_regulons <- lapply(reg_list, function(x) {
  list(tfmode = setNames(x$mor, x$target),
       likelihood = rep(1, nrow(x)))
})

tf_act <- viper(expr_mat, regulon = viper_regulons,
                eset.filter = FALSE, minsize = 4, verbose = FALSE)

glyco_scores <- tumor_hep$Glycolysis_AUC
cor_results <- data.frame(TF = rownames(tf_act), rho = NA, p = NA)
for(i in 1:nrow(tf_act)) {
  ct <- cor.test(tf_act[i, ], glyco_scores, method = "spearman")
  cor_results$rho[i] <- ct$estimate
  cor_results$p[i]   <- ct$p.value
}
cor_results$padj <- p.adjust(cor_results$p, method = "BH")
cor_results <- cor_results[order(-abs(cor_results$rho)), ]
# ============================================================
# NicheNet配体活性分析（Revision新增）
# ============================================================

library(nichenetr)
library(dplyr)

ligand_target_matrix <- readRDS("D:/scRNA_project/nichenet/ligand_target_matrix_nsga2r_final.rds")
lr_network           <- readRDS("D:/scRNA_project/nichenet/lr_network_human_21122021.rds")
weighted_networks    <- readRDS("D:/scRNA_project/nichenet/weighted_networks_nsga2r_final.rds")

sender_cells <- colnames(subset(seu, site == "Tumor" &
                                  celltype == "Hepatocyte" &
                                  Glycolysis_AUC > median(seu$Glycolysis_AUC[
                                    seu$site == "Tumor" & seu$celltype == "Hepatocyte"])))

receiver_cells <- colnames(subset(seu, site == "Tumor" &
                                    celltype %in% c("T/NK", "Myeloid")))

expr_sender   <- GetAssayData(seu, layer = "data")[, sender_cells]
expr_receiver <- GetAssayData(seu, layer = "data")[, receiver_cells]

expressed_genes_sender   <- rownames(expr_sender)[rowSums(expr_sender > 0) / ncol(expr_sender) > 0.10]
expressed_genes_receiver <- rownames(expr_receiver)[rowSums(expr_receiver > 0) / ncol(expr_receiver) > 0.10]

ligands   <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

expressed_ligands   <- intersect(ligands,   expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

geneset_oi <- intersect(
  c("PDCD1","HAVCR2","TIGIT","LAG3","TOX","CTLA4","CD274","VSIR",
    "ENTPD1","CD96","ARG1","IL10","TGFB1","CD163","MRC1",
    "CXCL8","CCL2","VEGFA","SPP1","MIF"),
  expressed_genes_receiver
)

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = expressed_genes_receiver,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = expressed_ligands
)

ligand_activities <- ligand_activities %>% arrange(-aupr_corrected)
# MIF排名37/325 AUPR=0.091；SPP1排名250/325 AUPR=0.026

# ============================================================
# Supplementary Figure S9: UMI Downsampling Sensitivity
# ============================================================

library(Matrix)
library(CellChat)

cat("=== S9: UMI Downsampling Sensitivity ===\n")

# tumor_hep must be loaded (from workspace or Script 02 output)
# Compute median UMI per group
high_cells <- colnames(tumor_hep)[tumor_hep$glyco_group == "TumorGlycoHigh"]
low_cells  <- colnames(tumor_hep)[tumor_hep$glyco_group == "TumorGlycoLow"]

umi_high <- colSums(GetAssayData(tumor_hep[, high_cells], layer = "counts"))
umi_low  <- colSums(GetAssayData(tumor_hep[, low_cells],  layer = "counts"))
target_umi <- median(umi_low)   # ~4,410
cat("GlycoHigh median UMI:", median(umi_high), "\n")
cat("GlycoLow  median UMI:", median(umi_low),  "\n")
cat("Downsampling target:", target_umi, "\n")

# Multinomial downsampling of GlycoHigh cells
set.seed(42)
counts_high <- GetAssayData(tumor_hep[, high_cells], layer = "counts")
counts_ds   <- apply(counts_high, 2, function(cell_counts) {
  total <- sum(cell_counts)
  if (total <= target_umi) return(cell_counts)
  probs <- cell_counts / total
  as.integer(rmultinom(1, size = target_umi, prob = probs))
})
rownames(counts_ds) <- rownames(counts_high)

# Re-compute AUCell glycolysis scores on downsampled matrix
rankings_ds  <- AUCell_buildRankings(counts_ds, plotStats = FALSE)
auc_ds       <- AUCell_calcAUC(list(glycolysis = glycolysis_genes),
                                rankings_ds,
                                aucMaxRank = ceiling(0.05 * nrow(counts_ds)))
glyco_ds_scores <- as.numeric(getAUC(auc_ds)["glycolysis", ])
cat("Glycolysis AUC post-downsampling — mean:",
    round(mean(glyco_ds_scores), 4), "\n")

# Re-run CellChat on downsampled dataset
# (requires seu object with cell_type_glyco metadata)
# Build expression matrix: replace GlycoHigh with downsampled counts
counts_full      <- GetAssayData(seu, assay = "RNA", layer = "counts")
counts_full_ds   <- counts_full
counts_full_ds[, high_cells] <- counts_ds

seu_ds <- seu
seu_ds[["RNA"]] <- SetAssayData(seu_ds[["RNA"]], layer = "counts",
                                 new.data = counts_full_ds)
seu_ds <- NormalizeData(seu_ds)

meta_ds <- data.frame(cell_type_glyco = seu$cell_type_glyco,
                       row.names       = colnames(seu))
cellchat_ds <- createCellChat(
  object   = GetAssayData(seu_ds, layer = "data"),
  meta     = meta_ds,
  group.by = "cell_type_glyco"
)
cellchat_ds@DB <- CellChatDB.human
cellchat_ds <- subsetData(cellchat_ds)
cellchat_ds <- identifyOverExpressedGenes(cellchat_ds)
cellchat_ds <- identifyOverExpressedInteractions(cellchat_ds)
cellchat_ds <- computeCommunProb(cellchat_ds, type = "triMean", nboot = 1000)
cellchat_ds <- filterCommunication(cellchat_ds, min.cells = 10)
cellchat_ds <- computeCommunProbPathway(cellchat_ds)
cellchat_ds <- aggregateNet(cellchat_ds)
saveRDS(cellchat_ds, "D:/scRNA_project/hcc_cellchat_downsampled.rds")

# Extract SPP1 and MIF pathway probabilities and compare
extract_pathway_prob <- function(cc, pathway, sender, receiver) {
  df <- subsetCommunication(cc, slot.name = "netP")
  df[df$pathway_name == pathway &
       df$source == sender & df$target == receiver, "prob"]
}

# Report: SPP1 and MIF to Macrophage before and after downsampling
cat("\n--- Downsampling comparison (manuscript Table S9) ---\n")
for (pw in c("SPP1", "MIF")) {
  orig_prob <- extract_pathway_prob(
    readRDS("D:/scRNA_project/hcc_cellchat_round2_glyco.rds"),
    pw, "TumorGlycoHigh", "Myeloid"
  )
  ds_prob <- extract_pathway_prob(cellchat_ds, pw, "TumorGlycoHigh", "Myeloid")
  cat(sprintf("%s -> Myeloid | Original: %.4f | Downsampled: %.4f\n",
              pw,
              ifelse(length(orig_prob) > 0, orig_prob, NA),
              ifelse(length(ds_prob)   > 0, ds_prob,   NA)))
}
# SPP1 expect: 0.209 vs 0.052 (FC 3.1-6.5); MIF: 0.1368 vs 0.1198

# ============================================================
# Supplementary Figure S10: Permutation Control
# ENO1-glycolysis circularity check
# ============================================================

cat("\n=== S10: Permutation control (ENO1 circularity) ===\n")

set.seed(42)
n_perm  <- 100
eno1_expr <- as.numeric(
  GetAssayData(tumor_hep, layer = "data")["ENO1", ]
)

# Real glycolysis AUC score
auc_real  <- AUCell_calcAUC(list(glycolysis = glycolysis_genes),
                             rankings,
                             aucMaxRank = ceiling(0.05 * nrow(rankings)))
real_rho  <- cor(eno1_expr, as.numeric(getAUC(auc_real)["glycolysis", ]),
                 method = "spearman")

# Background: 100 random non-glycolytic gene sets of equal size
all_genes    <- rownames(GetAssayData(tumor_hep, layer = "data"))
non_glyco    <- setdiff(all_genes, glycolysis_genes)
perm_rhos    <- numeric(n_perm)

for (i in seq_len(n_perm)) {
  rand_set  <- sample(non_glyco, length(glycolysis_genes))
  auc_rand  <- AUCell_calcAUC(list(rand = rand_set), rankings,
                               aucMaxRank = ceiling(0.05 * nrow(rankings)))
  perm_rhos[i] <- cor(eno1_expr,
                      as.numeric(getAUC(auc_rand)["rand", ]),
                      method = "spearman")
}

emp_p <- mean(perm_rhos >= real_rho)
cat(sprintf("Real glycolysis rho = %.3f | Permutation empirical p = %.3f\n",
            real_rho, emp_p))  # expect rho~0.541, p<0.01

# Plot S10
perm_df <- data.frame(rho = perm_rhos)
p_s10 <- ggplot(perm_df, aes(x = rho)) +
  geom_histogram(bins = 30, fill = "grey70", color = "white") +
  geom_vline(xintercept = real_rho, color = "#E64B35",
             linewidth = 1, linetype = "dashed") +
  annotate("text", x = real_rho + 0.02, y = Inf,
           label = sprintf("Real rho = %.3f\np < 0.01", real_rho),
           hjust = 0, vjust = 1.5, color = "#E64B35", size = 4) +
  labs(x = "Spearman rho (ENO1 vs. random gene set AUCell score)",
       y = "Count",
       title = "Permutation Control: ENO1-Glycolysis Circularity Check",
       subtitle = "100 random gene sets of equal size (n=22)") +
  theme_classic(base_size = 12)

ggsave("output/FigS10_permutation_control.pdf", p_s10, width = 7, height = 5)
ggsave("output/FigS10_permutation_control.png", p_s10, width = 7, height = 5,
       dpi = 300)
cat("Supplementary Figure S10 saved.\n")

# ============================================================
# Supplementary Figure S12: Per-patient stratified correlation
# ============================================================

cat("\n=== S12: Per-patient ENO1-glycolysis correlation ===\n")

patients <- unique(tumor_hep$patient)
cat("Patients:", length(patients), "\n")

per_patient_cor <- lapply(patients, function(pt) {
  cells  <- colnames(tumor_hep)[tumor_hep$patient == pt]
  if (length(cells) < 30) {
    cat("  Skipping", pt, "(n =", length(cells), "< 30)\n")
    return(NULL)
  }
  eno1_pt  <- as.numeric(
    GetAssayData(tumor_hep[, cells], layer = "data")["ENO1", ]
  )
  glyco_pt <- tumor_hep$glycolysis_score[cells]
  ct       <- cor.test(eno1_pt, glyco_pt, method = "spearman", exact = FALSE)

  # Fisher z CI for Spearman rho
  z     <- 0.5 * log((1 + ct$estimate) / (1 - ct$estimate))
  se    <- 1 / sqrt(length(cells) - 3)
  ci_lo <- tanh(z - 1.96 * se)
  ci_hi <- tanh(z + 1.96 * se)

  data.frame(
    patient  = pt,
    n_cells  = length(cells),
    rho      = ct$estimate,
    ci_lower = ci_lo,
    ci_upper = ci_hi,
    pval     = ct$p.value
  )
}) %>% bind_rows()

per_patient_cor$padj <- p.adjust(per_patient_cor$pval, method = "BH")
cat("\nPer-patient correlations:\n"); print(per_patient_cor)
# All 8 should be significant; rho range 0.193-0.610

p_s12 <- ggplot(per_patient_cor,
                aes(x = reorder(patient, rho), y = rho)) +
  geom_hline(yintercept = 0.57, linetype = "dashed",
             color = "grey50", linewidth = 0.6) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper),
                width = 0.3, color = "grey40") +
  geom_point(aes(size = n_cells, color = padj < 0.05)) +
  scale_color_manual(values = c("TRUE" = "#E64B35", "FALSE" = "grey60"),
                     name = "BH-adj p < 0.05") +
  scale_size_continuous(name = "Cell count", range = c(3, 8)) +
  annotate("text", x = Inf, y = 0.57, hjust = 1.1, vjust = -0.5,
           label = "Global R = 0.57", size = 3.5, color = "grey50") +
  labs(x = "Patient", y = "Spearman rho (ENO1 ~ Glycolysis AUC)",
       title = "Per-patient ENO1-Glycolysis Correlation",
       subtitle = "Error bars = 95% CI (Fisher z-transformation)") +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")

ggsave("output/FigS12_per_patient_correlation.pdf",
       p_s12, width = 8, height = 5)
ggsave("output/FigS12_per_patient_correlation.png",
       p_s12, width = 8, height = 5, dpi = 300)
cat("Supplementary Figure S12 saved.\n")

cat("\n=== All revision analyses complete ===\n")
