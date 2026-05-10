# ============================================================
# Script 07: LASSO Glycolysis Risk Score
# Section 2.10, Results Section 3.9
# Generates: Figure 7A-B, Supplementary Figure S14
# ============================================================

library(glmnet)
library(survival)
library(survminer)
library(ggplot2)

setwd("D:/scRNA_project")

# ── Load data ─────────────────────────────────────────────
expr     <- read.table("TCGA_LIHC_expression.txt",
                       header = TRUE, row.names = 1,
                       sep = "\t", check.names = FALSE)
clin     <- read.table("TCGA_LIHC_clinical.txt",
                       header = TRUE, sep = "\t", check.names = FALSE)

# ── Glycolysis gene set (22 genes) ────────────────────────
glyco_genes <- c("HK1","HK2","GPI","PFKL","PFKP","PFKM",
                 "ALDOA","ALDOB","ALDOC","TPI1","GAPDH","PGK1",
                 "PGAM1","ENO1","ENO2","PKM","LDHA","LDHB",
                 "SLC2A1","SLC2A3","PFKFB3","GCK")
glyco_genes <- glyco_genes[glyco_genes %in% rownames(expr)]
cat("Glycolysis genes available:", length(glyco_genes), "\n")

# Match samples with survival data
common_samples <- intersect(colnames(expr), clin$sample)
expr_sub       <- t(expr[glyco_genes, common_samples])
clin_sub       <- clin[match(common_samples, clin$sample), ]
keep           <- !is.na(clin_sub$OS.time) & !is.na(clin_sub$OS)
expr_sub       <- expr_sub[keep, ]
clin_sub       <- clin_sub[keep, ]
cat("Final sample count:", nrow(expr_sub), "\n")   # expect ~418

# ── LASSO Cox regression ───────────────────────────────────
set.seed(42)
cv_fit <- cv.glmnet(
  x      = as.matrix(expr_sub),
  y      = Surv(clin_sub$OS.time, clin_sub$OS),
  family = "cox",
  alpha  = 1,
  nfolds = 10
)

# Supplementary Figure S14
png("FigS14_LASSO_CV.png", width = 7, height = 5,
    units = "in", res = 300, bg = "white")
plot(cv_fit); title("10-fold CV: LASSO Cox")
dev.off()

# Extract non-zero coefficients at lambda.min
lasso_coef <- coef(cv_fit, s = "lambda.min")
nonzero    <- lasso_coef[lasso_coef[, 1] != 0, , drop = FALSE]
cat("Genes with non-zero coefficients at lambda.min:\n")
print(nonzero)   # expect 9 genes

# Verify lambda.1se comparison (should collapse to 0 genes)
lasso_coef_1se <- coef(cv_fit, s = "lambda.1se")
nonzero_1se    <- lasso_coef_1se[lasso_coef_1se[, 1] != 0, , drop = FALSE]
cat("Genes at lambda.1se:", nrow(nonzero_1se), "(expect 0)\n")

# ── Compute risk score ────────────────────────────────────
coef_vec   <- as.numeric(nonzero)
names(coef_vec) <- rownames(nonzero)

risk_score         <- as.matrix(expr_sub[, names(coef_vec)]) %*% coef_vec
clin_sub$risk_score <- as.numeric(risk_score)
clin_sub$risk_group <- ifelse(
  clin_sub$risk_score >= median(clin_sub$risk_score), "High", "Low"
)
cat("High risk:", sum(clin_sub$risk_group == "High"),
    "| Low risk:", sum(clin_sub$risk_group == "Low"), "\n")

write.csv(clin_sub[, c("sample","risk_score","risk_group","OS.time","OS")],
          "TCGA_LIHC_risk_scores.csv", row.names = FALSE)

# ── Kaplan-Meier — Figure 7A ──────────────────────────────
fit    <- survfit(Surv(OS.time, OS) ~ risk_group, data = clin_sub)
lr     <- survdiff(Surv(OS.time, OS) ~ risk_group, data = clin_sub)
pval   <- 1 - pchisq(lr$chisq, df = 1)
cat("KM log-rank p =", formatC(pval, format = "e", digits = 2), "\n")  # expect 2.66e-12

p_km <- ggsurvplot(
  fit, data = clin_sub,
  pval = TRUE, risk.table = TRUE,
  palette = c("#E64B35","#4DBBD5"),
  legend.labs = c("High risk","Low risk"),
  xlab = "Time (days)", ylab = "Overall survival probability",
  title = paste0("Glycolysis Risk Score: TCGA-LIHC (n=", nrow(clin_sub), ")"),
  ggtheme = theme_classic(base_size = 13)
)
ggsave("Fig7A_LASSO_KM.png", print(p_km), width = 8, height = 7, dpi = 300)

# ── Multivariate Cox — Figure 7B ─────────────────────────
clin_cox <- clin_sub %>%
  mutate(stage_binary = ifelse(
    grepl("III|IV", ajcc_pathologic_stage), "III-IV", "I-II"
  )) %>%
  filter(!is.na(stage_binary), !is.na(age_at_index), !is.na(gender))

cox_multi <- coxph(
  Surv(OS.time, OS) ~ risk_score + age_at_index + gender + stage_binary,
  data = clin_cox
)
cat("\nMultivariate Cox (HR for risk_score):\n")
print(summary(cox_multi)$conf.int)
# expect HR~3.214 for risk_score

# ── FIX: Bootstrap C-index (n=1000) ──────────────────────
cat("\nBootstrap C-index (n=1000)...\n")
set.seed(42)
n        <- nrow(clin_cox)
boot_c   <- numeric(1000)

for (b in seq_len(1000)) {
  idx      <- sample(n, n, replace = TRUE)
  boot_dat <- clin_cox[idx, ]
  tryCatch({
    boot_cox  <- coxph(
      Surv(OS.time, OS) ~ risk_score + age_at_index + gender + stage_binary,
      data = boot_dat
    )
    boot_c[b] <- summary(boot_cox)$concordance[1]
  }, error = function(e) {
    boot_c[b] <<- NA
  })
}
boot_c  <- boot_c[!is.na(boot_c)]
c_point <- summary(cox_multi)$concordance[1]
c_lower <- quantile(boot_c, 0.025)
c_upper <- quantile(boot_c, 0.975)
cat(sprintf("C-index = %.3f (95%% CI: %.3f–%.3f, n=%d bootstrap)\n",
            c_point, c_lower, c_upper, length(boot_c)))
# expect 0.685 (0.637–0.732)

# Forest plot
png("Fig7B_Cox_forest.png", width = 8, height = 5,
    units = "in", res = 300, bg = "white")
ggforest(cox_multi, data = clin_cox,
         main = "Multivariate Cox: Glycolysis Risk Score")
dev.off()

cat("\nLASSO coefficients:\n"); print(coef_vec)
cat("\n=== Script 07 complete ===\n")
