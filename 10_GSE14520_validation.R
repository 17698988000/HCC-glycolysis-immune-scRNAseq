# ============================================================
# 10_GSE14520_validation.R
# Purpose: Independent validation of glycolysis risk score
#          in GSE14520 microarray cohort
# Corresponds to: Methods Section 2.13, Results Section 3.12
# ============================================================

library(GEOquery)
library(hgu133plus2.db)
library(AnnotationDbi)
library(survival)
library(survminer)
library(ggplot2)

# ── LASSO coefficients from TCGA-LIHC training ────────────
lasso_coef <- c(
  ENO1   = 0.18113254,
  LDHA   = 0.11009435,
  TPI1   = 0.06877734,
  SLC2A1 = 0.73819068,
  PFKFB3 = 0.61190122,
  GPI    = 0.10475223,
  PFKL   = 0.04659883,
  PFKM   = 0.40162803,
  ALDOA  = 0.45702991
)

# ── Download GSE14520 ──────────────────────────────────────
cat("Downloading GSE14520...\n")
gse <- getGEO("GSE14520", GSEMatrix = TRUE, AnnotGPL = TRUE)

# Use GPL3921 platform (Affymetrix HT HG-U133A, 445 samples)
expr_main  <- exprs(gse[[1]])
cat("Expression matrix:", dim(expr_main), "\n")

# Load supplementary clinical data
# Download from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE14520
# File: GSE14520_Extra_Supplement.txt.gz
clin <- read.table(
  gzfile("D:/scRNA_project/GSE14520_Extra_Supplement.txt.gz"),
  header = TRUE, sep = "\t", stringsAsFactors = FALSE
)

# ── Probe-to-gene annotation ───────────────────────────────
probe_gene <- AnnotationDbi::select(
  hgu133plus2.db,
  keys    = rownames(expr_main),
  columns = c("PROBEID","SYMBOL"),
  keytype = "PROBEID"
)
probe_gene <- probe_gene[!is.na(probe_gene$SYMBOL), ]

# ── Extract 9 risk genes (mean across probes) ─────────────
risk_genes <- names(lasso_coef)
gene_expr  <- sapply(risk_genes, function(g) {
  probes <- probe_gene$PROBEID[probe_gene$SYMBOL == g]
  probes <- probes[probes %in% rownames(expr_main)]
  if (length(probes) == 0) return(rep(NA, ncol(expr_main)))
  if (length(probes) == 1) return(expr_main[probes, ])
  colMeans(expr_main[probes, ])
})
cat("Gene coverage:\n"); print(colSums(!is.na(gene_expr)))

# ── Compute risk scores ────────────────────────────────────
risk_score <- as.matrix(gene_expr) %*% lasso_coef

# Match with tumor clinical data
tumor_clin <- clin[clin$Tissue.Type == "Tumor", ]
common_id  <- intersect(colnames(expr_main), tumor_clin$Affy_GSM)
cat("Matched tumor samples:", length(common_id), "\n")

risk_df <- data.frame(
  sample     = common_id,
  risk_score = risk_score[common_id, 1],
  OS_time    = tumor_clin$Survival.months[match(common_id, tumor_clin$Affy_GSM)],
  OS_status  = tumor_clin$Survival.status[match(common_id, tumor_clin$Affy_GSM)]
)
risk_df    <- risk_df[complete.cases(risk_df), ]
risk_df$risk_group <- ifelse(
  risk_df$risk_score >= median(risk_df$risk_score), "High", "Low"
)
cat("Final sample count:", nrow(risk_df), "\n")

# ── KM analysis ───────────────────────────────────────────
lr    <- survdiff(Surv(OS_time, OS_status) ~ risk_group, data = risk_df)
pval  <- 1 - pchisq(lr$chisq, df = 1)
cox   <- coxph(Surv(OS_time, OS_status) ~ risk_group, data = risk_df)

cat("log-rank p-value:", pval, "\n")
cat("Cox HR (High vs Low):\n")
print(summary(cox)$conf.int)

# ── Plot FigS17 ────────────────────────────────────────────
fit <- survfit(Surv(OS_time, OS_status) ~ risk_group, data = risk_df)
p_km <- ggsurvplot(
  fit, data = risk_df,
  pval = TRUE, pval.method = TRUE,
  risk.table = TRUE, risk.table.height = 0.28,
  palette = c("#E64B35","#4DBBD5"),
  legend.labs = c("High risk","Low risk"),
  xlab = "Time (months)", ylab = "Overall survival probability",
  title = paste0("Glycolysis Risk Score Validation: GSE14520 (n = ", nrow(risk_df), ")"),
  ggtheme = theme_classic(base_size = 13),
  surv.median.line = "hv"
)

png("D:/scRNA_project/FigS17_GSE14520_KM_v2.png",
    width = 8, height = 7, units = "in", res = 300, bg = "white")
print(p_km)
dev.off()

pdf("D:/scRNA_project/FigS17_GSE14520_KM_v2.pdf",
    width = 8, height = 7, bg = "white")
print(p_km)
dev.off()

cat("GSE14520 validation complete.\n")
