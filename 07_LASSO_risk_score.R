# ============================================================
# 07_LASSO_risk_score.R
# Purpose: Construct nine-gene glycolysis risk score using
#          LASSO-penalized Cox regression in TCGA-LIHC
# Corresponds to: Methods Section 2.10
# ============================================================

library(glmnet)
library(survival)
library(survminer)
library(ggplot2)

# ── Load TCGA-LIHC data ────────────────────────────────────
# Expression data (STAR FPKM, log2-transformed)
expr <- read.table("D:/scRNA_project/TCGA_LIHC_expression.txt",
                   header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
# Clinical/survival data
clin <- read.table("D:/scRNA_project/TCGA_LIHC_clinical.txt",
                   header = TRUE, sep = "\t", check.names = FALSE)

# ── Glycolysis gene set ────────────────────────────────────
glyco_genes <- c("HK1","HK2","GPI","PFKL","PFKP","PFKM","ALDOA","ALDOB","ALDOC",
                 "TPI1","GAPDH","PGK1","PGAM1","ENO1","ENO2","PKM","LDHA","LDHB",
                 "SLC2A1","SLC2A3","PFKFB3","GCK")

# Subset to available genes
glyco_genes <- glyco_genes[glyco_genes %in% rownames(expr)]
cat("Glycolysis genes available:", length(glyco_genes), "\n")

# Match samples with survival data
common_samples <- intersect(colnames(expr), clin$sample)
expr_sub <- t(expr[glyco_genes, common_samples])
clin_sub  <- clin[match(common_samples, clin$sample), ]
cat("Samples with survival data:", nrow(expr_sub), "\n")

# Remove samples with missing survival
keep <- !is.na(clin_sub$OS.time) & !is.na(clin_sub$OS)
expr_sub <- expr_sub[keep, ]
clin_sub  <- clin_sub[keep, ]
cat("Final sample count:", nrow(expr_sub), "\n")

# ── LASSO Cox regression ───────────────────────────────────
set.seed(42)
cv_fit <- cv.glmnet(
  x      = as.matrix(expr_sub),
  y      = Surv(clin_sub$OS.time, clin_sub$OS),
  family = "cox",
  alpha  = 1,
  nfolds = 10
)

# Save cross-validation plot (Supplementary Figure S14)
png("D:/scRNA_project/FigS14_LASSO_CV.png", width = 7, height = 5,
    units = "in", res = 300, bg = "white")
plot(cv_fit)
title("10-fold CV: LASSO Cox Regression")
dev.off()

# Extract non-zero coefficients at lambda.min
lasso_coef <- coef(cv_fit, s = "lambda.min")
nonzero    <- lasso_coef[lasso_coef[, 1] != 0, , drop = FALSE]
cat("Genes with non-zero coefficients:\n")
print(nonzero)

# ── Compute risk score ─────────────────────────────────────
coef_vec   <- as.numeric(nonzero)
names(coef_vec) <- rownames(nonzero)

risk_score <- as.matrix(expr_sub[, names(coef_vec)]) %*% coef_vec
clin_sub$risk_score  <- as.numeric(risk_score)
clin_sub$risk_group  <- ifelse(clin_sub$risk_score >= median(clin_sub$risk_score),
                                "High", "Low")

# Save risk score table
write.csv(clin_sub[, c("sample","risk_score","risk_group","OS.time","OS")],
          "D:/scRNA_project/TCGA_LIHC_risk_scores.csv", row.names = FALSE)

# ── Kaplan-Meier analysis ──────────────────────────────────
fit <- survfit(Surv(OS.time, OS) ~ risk_group, data = clin_sub)
lr  <- survdiff(Surv(OS.time, OS) ~ risk_group, data = clin_sub)
pval <- 1 - pchisq(lr$chisq, df = 1)
cat("KM log-rank p-value:", pval, "\n")

p_km <- ggsurvplot(
  fit, data = clin_sub,
  pval = TRUE, risk.table = TRUE,
  palette = c("#E64B35","#4DBBD5"),
  legend.labs = c("High risk","Low risk"),
  xlab = "Time (days)", ylab = "Overall survival probability",
  title = "Glycolysis Risk Score: TCGA-LIHC",
  ggtheme = theme_classic(base_size = 13)
)

png("D:/scRNA_project/Fig7A_LASSO_KM.png", width = 8, height = 7,
    units = "in", res = 300, bg = "white")
print(p_km)
dev.off()

# ── Multivariate Cox regression ────────────────────────────
clin_sub$stage_binary <- ifelse(
  clin_sub$ajcc_pathologic_stage %in% c("Stage III","Stage IV"), "III-IV", "I-II"
)
cox_multi <- coxph(
  Surv(OS.time, OS) ~ risk_score + age_at_index + gender + stage_binary,
  data = clin_sub
)
cat("\nMultivariate Cox results:\n")
print(summary(cox_multi)$conf.int)
cat("C-index:", summary(cox_multi)$concordance[1], "\n")

# Forest plot
png("D:/scRNA_project/Fig7B_Cox_forest.png", width = 8, height = 5,
    units = "in", res = 300, bg = "white")
ggforest(cox_multi, data = clin_sub,
         main = "Multivariate Cox: Glycolysis Risk Score")
dev.off()

cat("LASSO risk score analysis complete.\n")
cat("LASSO coefficients:\n")
print(coef_vec)
