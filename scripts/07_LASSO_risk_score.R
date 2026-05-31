#!/usr/bin/env Rscript

# ============================================================
# 07_LASSO_risk_score.R
#
# Purpose:
#   Reproduce the finalized TCGA-LIHC four-gene glycolysis score
#   construction and survival analysis.
#
# Manuscript-aligned convention:
#   Dataset: TCGA-LIHC primary tumor samples
#   Expression: log2(FPKM + 1)
#   Survival-matched LASSO cohort: n = 365
#   Candidate genes: 22 curated glycolysis genes
#
# Locked final model:
#   lambda.min retained 4 genes:
#     TPI1, ENO1, LDHA, SLC2A1
#
#   lambda.1se selected a null model.
#
#   Four-gene glycolysis score =
#     0.3041908 * TPI1 +
#     0.9639654 * ENO1 +
#     1.3404374 * LDHA +
#     0.2424239 * SLC2A1
#
# Outputs:
#   results/07_LASSO_risk_score/
#     07_TCGA_four_gene_score_source_data.csv
#     07_LASSO_lambda_min_coefficients.csv
#     07_LASSO_lambda_1se_coefficients.csv
#     07_LASSO_repeated_cv_selection_frequency.csv
#     07_TCGA_multivariable_cox_source.csv
#     07_TCGA_model_comparison_source.csv
#     07_LASSO_risk_score_QC_check.csv
#     FigS14_LASSO_CV.pdf/png
#     Fig6A_TCGA_four_gene_score_KM.pdf/png
#     Fig6B_TCGA_four_gene_score_multivariable_Cox.pdf/png
# ============================================================

# ============================================================
# Package checks
# ============================================================

required_packages <- c(
  "glmnet",
  "survival",
  "survminer",
  "ggplot2",
  "dplyr",
  "tibble",
  "readr",
  "stringr"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    "\nInstall them before running this script."
  )
}

suppressPackageStartupMessages({
  library(glmnet)
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(readr)
  library(stringr)
})

# ---------- avoid masking ----------
select      <- dplyr::select
mutate      <- dplyr::mutate
filter      <- dplyr::filter
arrange     <- dplyr::arrange
summarise   <- dplyr::summarise
group_by    <- dplyr::group_by
ungroup     <- dplyr::ungroup
left_join   <- dplyr::left_join
inner_join  <- dplyr::inner_join
distinct    <- dplyr::distinct
case_when   <- dplyr::case_when
n_distinct  <- dplyr::n_distinct

# ============================================================
# User-facing inputs
# ============================================================

expr_path <- "TCGA_LIHC_expression.txt"
clin_path <- "TCGA_LIHC_clinical.txt"

outdir <- "results/07_LASSO_risk_score"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Set TRUE only if TCGA_LIHC_expression.txt is already log2(FPKM + 1).
# Final manuscript convention is log2(FPKM + 1), so the default assumes
# the input file contains untransformed FPKM-like values.
expr_values_are_already_log2 <- FALSE

restrict_to_tcga_primary_tumor <- TRUE

run_repeated_cv <- TRUE
n_repeated_cv <- 100L

set.seed(42)

# ============================================================
# Locked manuscript constants
# ============================================================

glyco_genes_22 <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "PFKM",
  "ALDOA", "ALDOB", "ALDOC", "TPI1", "GAPDH", "PGK1",
  "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB",
  "SLC2A1", "SLC2A3", "PFKFB3", "GCK"
)

final_coef <- c(
  TPI1   = 0.3041908,
  ENO1   = 0.9639654,
  LDHA   = 1.3404374,
  SLC2A1 = 0.2424239
)

expected_lasso_n <- 365L
expected_multivariable_cox_n <- 341L
expected_high_n <- 183L
expected_low_n <- 182L
expected_lambda_min_genes <- names(final_coef)
expected_lambda_1se_n <- 0L

# ============================================================
# Helper functions
# ============================================================

standardize_barcode <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("\\.", "-", x)
  x
}

tcga_patient_barcode <- function(x) {
  x <- standardize_barcode(x)
  out <- ifelse(grepl("^TCGA-", x) & nchar(x) >= 12, substr(x, 1, 12), x)
  out
}

tcga_sample_barcode_16 <- function(x) {
  x <- standardize_barcode(x)
  out <- ifelse(grepl("^TCGA-", x) & nchar(x) >= 16, substr(x, 1, 16), x)
  out
}

tcga_sample_type_code <- function(x) {
  x <- standardize_barcode(x)
  out <- rep(NA_character_, length(x))
  idx <- grepl("^TCGA-", x) & nchar(x) >= 15
  out[idx] <- substr(x[idx], 14, 15)
  out
}

find_col <- function(df, candidates, required = TRUE, label = "column") {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) {
    return(hit[1])
  }

  if (required) {
    stop(
      "Could not find required ",
      label,
      ". Tried: ",
      paste(candidates, collapse = ", "),
      "\nAvailable columns:\n",
      paste(colnames(df), collapse = ", ")
    )
  }

  NULL
}

parse_os_event <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(as.integer(x))
  }

  x_chr <- tolower(trimws(as.character(x)))

  dplyr::case_when(
    x_chr %in% c("1", "event", "dead", "deceased", "death", "true") ~ 1L,
    x_chr %in% c("0", "censored", "alive", "living", "false") ~ 0L,
    grepl("dead", x_chr) ~ 1L,
    grepl("alive", x_chr) ~ 0L,
    TRUE ~ NA_integer_
  )
}

parse_stage <- function(x) {
  x_chr <- toupper(trimws(as.character(x)))
  x_chr <- gsub("AJCC", "", x_chr)
  x_chr <- gsub("PATHOLOGIC", "", x_chr)
  x_chr <- gsub("PATHOLOGICAL", "", x_chr)
  x_chr <- gsub("STAGE", "", x_chr)
  x_chr <- trimws(x_chr)

  dplyr::case_when(
    is.na(x_chr) | x_chr == "" ~ NA_character_,
    grepl("IV", x_chr) ~ "IV",
    grepl("III", x_chr) ~ "III",
    grepl("II", x_chr) ~ "II",
    grepl("I", x_chr) ~ "I",
    TRUE ~ NA_character_
  )
}

cox_to_table <- function(fit) {
  s <- summary(fit)

  tibble::tibble(
    term = rownames(s$coefficients),
    beta = s$coefficients[, "coef"],
    HR = s$conf.int[, "exp(coef)"],
    lower_95_CI = s$conf.int[, "lower .95"],
    upper_95_CI = s$conf.int[, "upper .95"],
    p_value = s$coefficients[, "Pr(>|z|)"]
  )
}

selected_genes_from_cv <- function(cv_fit, s_value) {
  coef_mat <- as.matrix(coef(cv_fit, s = s_value))
  genes <- rownames(coef_mat)[coef_mat[, 1] != 0]
  genes
}

coef_table_from_cv <- function(cv_fit, s_value) {
  coef_mat <- as.matrix(coef(cv_fit, s = s_value))
  tibble::tibble(
    gene = rownames(coef_mat),
    coefficient = as.numeric(coef_mat[, 1])
  ) %>%
    filter(coefficient != 0) %>%
    arrange(desc(abs(coefficient)))
}

# ============================================================
# Load expression and clinical data
# ============================================================

if (!file.exists(expr_path)) {
  stop("Expression file not found: ", expr_path)
}

if (!file.exists(clin_path)) {
  stop("Clinical file not found: ", clin_path)
}

expr_raw <- read.table(
  expr_path,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

clin <- read.table(
  clin_path,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

if (nrow(expr_raw) == 0 || ncol(expr_raw) == 0) {
  stop("Expression matrix is empty.")
}

if (nrow(clin) == 0) {
  stop("Clinical table is empty.")
}

expr_mat <- as.matrix(expr_raw)
storage.mode(expr_mat) <- "numeric"

if (anyNA(expr_mat)) {
  stop("NA values detected after converting expression matrix to numeric.")
}

if (any(expr_mat < 0, na.rm = TRUE)) {
  stop("Negative expression values detected. This script expects FPKM-like or log2(FPKM+1)-like non-negative values.")
}

if (expr_values_are_already_log2) {
  expr_log <- expr_mat
  expr_transform <- "input already log2(FPKM + 1)"
} else {
  expr_log <- log2(expr_mat + 1)
  expr_transform <- "log2(input + 1)"
}

message("Expression transform used: ", expr_transform)

# ============================================================
# Restrict expression matrix to TCGA primary solid tumor if possible
# ============================================================

sample_type_codes <- tcga_sample_type_code(colnames(expr_log))
primary_filter_applied <- FALSE

if (restrict_to_tcga_primary_tumor && any(!is.na(sample_type_codes))) {
  primary_cols <- colnames(expr_log)[sample_type_codes == "01"]

  if (length(primary_cols) == 0) {
    stop("TCGA sample-type codes were detected, but no primary solid tumor code '01' samples were found.")
  }

  expr_log <- expr_log[, primary_cols, drop = FALSE]
  primary_filter_applied <- TRUE
  message("Restricted expression matrix to TCGA primary solid tumor samples: ", ncol(expr_log))
} else {
  message("No TCGA sample-type codes detected in expression column names; assuming expression matrix is already tumor-focused.")
}

# ============================================================
# Clinical column resolution
# ============================================================

if (!any(c("sample", "Sample", "barcode", "Barcode", "submitter_id", "bcr_patient_barcode", "patient") %in% colnames(clin))) {
  clin <- tibble::rownames_to_column(as.data.frame(clin), "sample")
}

sample_col <- find_col(
  clin,
  c("sample", "Sample", "barcode", "Barcode", "submitter_id", "bcr_patient_barcode", "patient"),
  required = TRUE,
  label = "sample/barcode column"
)

os_time_col <- find_col(
  clin,
  c("OS.time", "OS_time", "OS.time.days", "OS_days", "overall_survival_time", "days_to_death_or_last_follow_up"),
  required = FALSE,
  label = "overall survival time column"
)

os_event_col <- find_col(
  clin,
  c("OS", "OS_event", "OS.status", "OS_status", "overall_survival_event", "vital_status"),
  required = FALSE,
  label = "overall survival event column"
)

age_col <- find_col(
  clin,
  c("age_at_index", "age", "age_at_diagnosis", "age_at_initial_pathologic_diagnosis"),
  required = FALSE,
  label = "age column"
)

sex_col <- find_col(
  clin,
  c("gender", "sex", "Gender", "Sex"),
  required = FALSE,
  label = "sex/gender column"
)

stage_col <- find_col(
  clin,
  c("ajcc_pathologic_stage", "pathologic_stage", "pathologic.stage", "stage", "tumor_stage"),
  required = FALSE,
  label = "AJCC stage column"
)

# Survival time fallback using TCGA-style fields if OS.time is absent.
if (is.null(os_time_col)) {
  death_col <- find_col(
    clin,
    c("days_to_death", "days.to.death"),
    required = FALSE,
    label = "days_to_death column"
  )

  follow_col <- find_col(
    clin,
    c("days_to_last_follow_up", "days_to_last_followup", "days.to.last.follow.up"),
    required = FALSE,
    label = "days_to_last_follow_up column"
  )

  if (is.null(death_col) || is.null(follow_col)) {
    stop("Could not derive OS time. Need OS.time or both days_to_death and days_to_last_follow_up.")
  }

  clin$OS_time_derived <- ifelse(
    !is.na(suppressWarnings(as.numeric(clin[[death_col]]))),
    suppressWarnings(as.numeric(clin[[death_col]])),
    suppressWarnings(as.numeric(clin[[follow_col]]))
  )

  os_time_col <- "OS_time_derived"
}

if (is.null(os_event_col)) {
  vital_col <- find_col(
    clin,
    c("vital_status", "Vital.status", "vital.status"),
    required = TRUE,
    label = "vital status column"
  )
  clin$OS_event_derived <- parse_os_event(clin[[vital_col]])
  os_event_col <- "OS_event_derived"
}

# Age fallback from days_to_birth if explicit age is absent.
if (is.null(age_col)) {
  birth_col <- find_col(
    clin,
    c("days_to_birth", "days.to.birth"),
    required = FALSE,
    label = "days_to_birth column"
  )

  if (!is.null(birth_col)) {
    clin$age_derived <- abs(suppressWarnings(as.numeric(clin[[birth_col]]))) / 365.25
    age_col <- "age_derived"
  }
}

if (is.null(age_col)) {
  warning("Age column not found. Multivariable Cox model will fail unless this is corrected.")
}

if (is.null(sex_col)) {
  warning("Sex/gender column not found. Multivariable Cox model will fail unless this is corrected.")
}

if (is.null(stage_col)) {
  warning("AJCC stage column not found. Multivariable Cox model will fail unless this is corrected.")
}

# ============================================================
# Match expression samples to clinical records
# ============================================================

expr_keys <- tibble::tibble(
  expr_col = colnames(expr_log),
  exact_key = standardize_barcode(colnames(expr_log)),
  sample16_key = tcga_sample_barcode_16(colnames(expr_log)),
  patient12_key = tcga_patient_barcode(colnames(expr_log))
)

clin_keys <- clin %>%
  mutate(
    clin_row_id = dplyr::row_number(),
    clinical_sample = .data[[sample_col]],
    exact_key = standardize_barcode(.data[[sample_col]]),
    sample16_key = tcga_sample_barcode_16(.data[[sample_col]]),
    patient12_key = tcga_patient_barcode(.data[[sample_col]])
  )

match_stats <- tibble::tibble(
  key_type = c("exact_key", "sample16_key", "patient12_key"),
  n_matches = c(
    length(intersect(expr_keys$exact_key, clin_keys$exact_key)),
    length(intersect(expr_keys$sample16_key, clin_keys$sample16_key)),
    length(intersect(expr_keys$patient12_key, clin_keys$patient12_key))
  )
)

match_key_type <- match_stats$key_type[which.max(match_stats$n_matches)]
message("Clinical matching key used: ", match_key_type)
print(match_stats)

if (max(match_stats$n_matches) == 0) {
  stop("No overlap between expression sample IDs and clinical sample IDs.")
}

expr_match <- expr_keys %>%
  mutate(match_key = .data[[match_key_type]]) %>%
  filter(match_key %in% clin_keys[[match_key_type]]) %>%
  arrange(expr_col) %>%
  group_by(match_key) %>%
  slice(1) %>%
  ungroup()

clin_match <- clin_keys %>%
  mutate(match_key = .data[[match_key_type]]) %>%
  filter(match_key %in% expr_match$match_key) %>%
  arrange(clin_row_id) %>%
  group_by(match_key) %>%
  slice(1) %>%
  ungroup()

matched <- inner_join(
  expr_match,
  clin_match %>% select(match_key, clin_row_id, clinical_sample),
  by = "match_key"
) %>%
  arrange(expr_col)

clin_matched <- clin[matched$clin_row_id, , drop = FALSE]
clin_matched$expr_col <- matched$expr_col
clin_matched$match_key <- matched$match_key
clin_matched$clinical_sample <- matched$clinical_sample

# ============================================================
# Prepare survival and clinical covariates
# ============================================================

surv_df <- clin_matched %>%
  mutate(
    OS_time = suppressWarnings(as.numeric(.data[[os_time_col]])),
    OS_event = parse_os_event(.data[[os_event_col]])
  )

if (!is.null(age_col)) {
  surv_df$age <- suppressWarnings(as.numeric(surv_df[[age_col]]))
}

if (!is.null(sex_col)) {
  surv_df$sex <- factor(tolower(trimws(as.character(surv_df[[sex_col]]))))
}

if (!is.null(stage_col)) {
  surv_df$AJCC_stage <- factor(
    parse_stage(surv_df[[stage_col]]),
    levels = c("I", "II", "III", "IV")
  )
}

keep_survival <- !is.na(surv_df$OS_time) &
  !is.na(surv_df$OS_event) &
  surv_df$OS_time > 0

surv_df <- surv_df[keep_survival, , drop = FALSE]

message("Survival-matched samples after filtering: ", nrow(surv_df))

# ============================================================
# Prepare expression matrix for LASSO and final score
# ============================================================

missing_glyco_genes <- setdiff(glyco_genes_22, rownames(expr_log))
missing_final_genes <- setdiff(names(final_coef), rownames(expr_log))

if (length(missing_glyco_genes) > 0) {
  stop(
    "Missing 22-gene glycolysis candidate gene(s): ",
    paste(missing_glyco_genes, collapse = ", ")
  )
}

if (length(missing_final_genes) > 0) {
  stop(
    "Missing final four-gene score gene(s): ",
    paste(missing_final_genes, collapse = ", ")
  )
}

expr_lasso <- t(expr_log[glyco_genes_22, surv_df$expr_col, drop = FALSE])
expr_final <- t(expr_log[names(final_coef), surv_df$expr_col, drop = FALSE])

if (anyNA(expr_lasso) || anyNA(expr_final)) {
  stop("NA detected in model expression matrix.")
}

# ============================================================
# LASSO Cox model for Supplementary Figure S14
# ============================================================

set.seed(42)

cv_fit <- cv.glmnet(
  x = as.matrix(expr_lasso),
  y = survival::Surv(surv_df$OS_time, surv_df$OS_event),
  family = "cox",
  alpha = 1,
  nfolds = 10,
  standardize = TRUE
)

lambda_min_coef <- coef_table_from_cv(cv_fit, "lambda.min")
lambda_1se_coef <- coef_table_from_cv(cv_fit, "lambda.1se")

lambda_min_genes <- lambda_min_coef$gene
lambda_1se_genes <- lambda_1se_coef$gene

write_csv(
  lambda_min_coef,
  file.path(outdir, "07_LASSO_lambda_min_coefficients.csv")
)

write_csv(
  lambda_1se_coef,
  file.path(outdir, "07_LASSO_lambda_1se_coefficients.csv")
)

# ============================================================
# Repeated 10-fold CV selection stability
# ============================================================

selection_frequency <- tibble::tibble(
  gene = glyco_genes_22,
  lambda_min_selected_n = 0L,
  lambda_1se_selected_n = 0L
)

if (run_repeated_cv) {
  message("Running repeated 10-fold CV stability analysis: ", n_repeated_cv, " repetitions")

  repeated_records <- list()

  for (i in seq_len(n_repeated_cv)) {
    set.seed(1000 + i)

    cv_i <- cv.glmnet(
      x = as.matrix(expr_lasso),
      y = survival::Surv(surv_df$OS_time, surv_df$OS_event),
      family = "cox",
      alpha = 1,
      nfolds = 10,
      standardize = TRUE
    )

    sel_min_i <- selected_genes_from_cv(cv_i, "lambda.min")
    sel_1se_i <- selected_genes_from_cv(cv_i, "lambda.1se")

    repeated_records[[i]] <- tibble::tibble(
      repeat_id = i,
      lambda = c(rep("lambda.min", length(sel_min_i)), rep("lambda.1se", length(sel_1se_i))),
      gene = c(sel_min_i, sel_1se_i)
    )
  }

  repeated_selection <- dplyr::bind_rows(repeated_records)

  selection_frequency <- selection_frequency %>%
    left_join(
      repeated_selection %>%
        filter(lambda == "lambda.min") %>%
        count(gene, name = "lambda_min_selected_n"),
      by = "gene",
      suffix = c("", ".new")
    ) %>%
    left_join(
      repeated_selection %>%
        filter(lambda == "lambda.1se") %>%
        count(gene, name = "lambda_1se_selected_n"),
      by = "gene",
      suffix = c("", ".new")
    ) %>%
    mutate(
      lambda_min_selected_n = dplyr::coalesce(lambda_min_selected_n.new, 0L),
      lambda_1se_selected_n = dplyr::coalesce(lambda_1se_selected_n.new, 0L)
    ) %>%
    select(gene, lambda_min_selected_n, lambda_1se_selected_n) %>%
    arrange(desc(lambda_min_selected_n), gene)

  write_csv(
    repeated_selection,
    file.path(outdir, "07_LASSO_repeated_cv_raw_selection.csv")
  )
}

write_csv(
  selection_frequency,
  file.path(outdir, "07_LASSO_repeated_cv_selection_frequency.csv")
)

# ============================================================
# Compute locked final four-gene score
# ============================================================

four_gene_score <- as.numeric(as.matrix(expr_final[, names(final_coef), drop = FALSE]) %*% final_coef)
score_median <- median(four_gene_score, na.rm = TRUE)

score_df <- surv_df %>%
  mutate(
    TPI1 = expr_final[, "TPI1"],
    ENO1 = expr_final[, "ENO1"],
    LDHA = expr_final[, "LDHA"],
    SLC2A1 = expr_final[, "SLC2A1"],
    four_gene_score = four_gene_score,
    four_gene_score_z = as.numeric(scale(four_gene_score)),
    score_group = ifelse(four_gene_score >= score_median, "High", "Low"),
    score_group = factor(score_group, levels = c("Low", "High"))
  )

score_group_counts <- table(score_df$score_group)
low_n <- as.integer(score_group_counts[["Low"]])
high_n <- as.integer(score_group_counts[["High"]])

message("Four-gene score group counts: Low = ", low_n, "; High = ", high_n)

source_data <- score_df %>%
  select(
    expr_col,
    match_key,
    clinical_sample,
    OS_time,
    OS_event,
    age,
    sex,
    AJCC_stage,
    TPI1,
    ENO1,
    LDHA,
    SLC2A1,
    four_gene_score,
    four_gene_score_z,
    score_group
  )

write_csv(
  source_data,
  file.path(outdir, "07_TCGA_four_gene_score_source_data.csv")
)

# ============================================================
# Kaplan-Meier analysis
# ============================================================

km_fit <- survival::survfit(
  survival::Surv(OS_time, OS_event) ~ score_group,
  data = score_df
)

km_diff <- survival::survdiff(
  survival::Surv(OS_time, OS_event) ~ score_group,
  data = score_df
)

km_p <- 1 - pchisq(km_diff$chisq, df = 1)

km_source <- tibble::tibble(
  n = nrow(score_df),
  low_n = low_n,
  high_n = high_n,
  median_score_cutoff = score_median,
  logrank_chisq = km_diff$chisq,
  logrank_p = km_p
)

write_csv(
  km_source,
  file.path(outdir, "07_TCGA_KM_source.csv")
)

# ============================================================
# Multivariable Cox model
# ============================================================

cox_df <- score_df %>%
  filter(
    !is.na(OS_time),
    !is.na(OS_event),
    !is.na(four_gene_score_z),
    !is.na(age),
    !is.na(sex),
    !is.na(AJCC_stage)
  ) %>%
  mutate(
    sex = droplevels(sex),
    AJCC_stage = droplevels(AJCC_stage)
  )

message("Multivariable Cox complete-case n: ", nrow(cox_df))

cox_clinical <- survival::coxph(
  survival::Surv(OS_time, OS_event) ~ age + sex + AJCC_stage,
  data = cox_df,
  x = TRUE,
  y = TRUE
)

cox_plus_score <- survival::coxph(
  survival::Surv(OS_time, OS_event) ~ four_gene_score_z + age + sex + AJCC_stage,
  data = cox_df,
  x = TRUE,
  y = TRUE
)

cox_source <- cox_to_table(cox_plus_score)

write_csv(
  cox_source,
  file.path(outdir, "07_TCGA_multivariable_cox_source.csv")
)

lrt <- anova(cox_clinical, cox_plus_score, test = "LRT")

model_comparison <- tibble::tibble(
  model = c("clinical_only", "clinical_plus_four_gene_score"),
  n = c(nrow(cox_df), nrow(cox_df)),
  AIC = c(AIC(cox_clinical), AIC(cox_plus_score)),
  C_index = c(
    summary(cox_clinical)$concordance[1],
    summary(cox_plus_score)$concordance[1]
  ),
  LRT_chisq_vs_previous = c(NA_real_, lrt$Chisq[2]),
  LRT_p_vs_previous = c(NA_real_, lrt$`Pr(>|Chi|)`[2])
)

write_csv(
  model_comparison,
  file.path(outdir, "07_TCGA_model_comparison_source.csv")
)

# ============================================================
# QC before plotting
# ============================================================

lambda_min_genes_match <- identical(
  sort(lambda_min_genes),
  sort(expected_lambda_min_genes)
)

lambda_1se_is_null <- length(lambda_1se_genes) == expected_lambda_1se_n

expected_genes_selected_100 <- TRUE
lambda_1se_null_in_repeated_cv <- TRUE

if (run_repeated_cv) {
  expected_genes_selected_100 <- selection_frequency %>%
    filter(gene %in% expected_lambda_min_genes) %>%
    summarise(ok = all(lambda_min_selected_n == n_repeated_cv)) %>%
    pull(ok)

  lambda_1se_null_in_repeated_cv <- all(selection_frequency$lambda_1se_selected_n == 0)
}

qc_check <- tibble::tibble(
  check_name = c(
    "expression_file_exists",
    "clinical_file_exists",
    "expression_transform",
    "primary_tumor_filter_applied_or_assumed",
    "all_22_glycolysis_genes_present",
    "all_final_4_genes_present",
    "TCGA_survival_matched_n",
    "lambda.min_selected_genes",
    "lambda.1se_selected_gene_n",
    "locked_score_formula_used",
    "score_group_Low_n",
    "score_group_High_n",
    "multivariable_Cox_complete_case_n",
    "repeated_CV_expected_4_genes_selected_100_of_100",
    "repeated_CV_lambda_1se_null_all_repeats",
    "overall_status"
  ),
  expected = c(
    "TRUE",
    "TRUE",
    "log2(FPKM + 1)",
    "TRUE or already tumor-focused",
    paste(glyco_genes_22, collapse = ";"),
    paste(names(final_coef), collapse = ";"),
    as.character(expected_lasso_n),
    paste(expected_lambda_min_genes, collapse = ";"),
    as.character(expected_lambda_1se_n),
    "0.3041908*TPI1 + 0.9639654*ENO1 + 1.3404374*LDHA + 0.2424239*SLC2A1",
    as.character(expected_low_n),
    as.character(expected_high_n),
    as.character(expected_multivariable_cox_n),
    ifelse(run_repeated_cv, as.character(n_repeated_cv), "not run"),
    ifelse(run_repeated_cv, "TRUE", "not run"),
    "PASS"
  ),
  observed = c(
    as.character(file.exists(expr_path)),
    as.character(file.exists(clin_path)),
    expr_transform,
    as.character(primary_filter_applied || !any(!is.na(sample_type_codes))),
    paste(setdiff(glyco_genes_22, missing_glyco_genes), collapse = ";"),
    paste(setdiff(names(final_coef), missing_final_genes), collapse = ";"),
    as.character(nrow(score_df)),
    paste(lambda_min_genes, collapse = ";"),
    as.character(length(lambda_1se_genes)),
    paste(
      paste0(names(final_coef), "=", final_coef),
      collapse = ";"
    ),
    as.character(low_n),
    as.character(high_n),
    as.character(nrow(cox_df)),
    ifelse(run_repeated_cv, as.character(expected_genes_selected_100), "not run"),
    ifelse(run_repeated_cv, as.character(lambda_1se_null_in_repeated_cv), "not run"),
    "PASS"
  ),
  pass = c(
    file.exists(expr_path),
    file.exists(clin_path),
    TRUE,
    primary_filter_applied || !any(!is.na(sample_type_codes)),
    length(missing_glyco_genes) == 0,
    length(missing_final_genes) == 0,
    nrow(score_df) == expected_lasso_n,
    lambda_min_genes_match,
    lambda_1se_is_null,
    TRUE,
    low_n == expected_low_n,
    high_n == expected_high_n,
    nrow(cox_df) == expected_multivariable_cox_n,
    ifelse(run_repeated_cv, expected_genes_selected_100, TRUE),
    ifelse(run_repeated_cv, lambda_1se_null_in_repeated_cv, TRUE),
    TRUE
  )
)

qc_check$pass[qc_check$check_name == "overall_status"] <- all(qc_check$pass[-nrow(qc_check)])
qc_check$observed[qc_check$check_name == "overall_status"] <- ifelse(
  all(qc_check$pass[-nrow(qc_check)]),
  "PASS",
  "FAIL"
)

write_csv(
  qc_check,
  file.path(outdir, "07_LASSO_risk_score_QC_check.csv")
)

if (!all(qc_check$pass)) {
  print(qc_check)
  stop(
    "QC failed. Source CSV and QC CSV were written, but figures were not generated. ",
    "Check sample matching, expression transformation, TCGA primary-tumor filtering, and LASSO preprocessing."
  )
}

message("QC PASS. Generating figures.")

# ============================================================
# Figure S14: LASSO CV curve
# ============================================================

pdf(
  file.path(outdir, "FigS14_LASSO_CV.pdf"),
  width = 6,
  height = 5,
  useDingbats = FALSE
)
plot(cv_fit)
dev.off()

png(
  file.path(outdir, "FigS14_LASSO_CV.png"),
  width = 6,
  height = 5,
  units = "in",
  res = 300,
  bg = "white"
)
plot(cv_fit)
dev.off()

# ============================================================
# Figure 6A: Kaplan-Meier plot
# ============================================================

km_plot <- survminer::ggsurvplot(
  km_fit,
  data = score_df,
  pval = TRUE,
  risk.table = TRUE,
  conf.int = FALSE,
  legend.title = "Four-gene score",
  legend.labs = c("Low", "High"),
  xlab = "Time (days)",
  ylab = "Overall survival probability",
  title = NULL,
  ggtheme = ggplot2::theme_classic(base_size = 13)
)

pdf(
  file.path(outdir, "Fig6A_TCGA_four_gene_score_KM.pdf"),
  width = 7,
  height = 6,
  useDingbats = FALSE
)
print(km_plot)
dev.off()

png(
  file.path(outdir, "Fig6A_TCGA_four_gene_score_KM.png"),
  width = 7,
  height = 6,
  units = "in",
  res = 300,
  bg = "white"
)
print(km_plot)
dev.off()

# ============================================================
# Figure 6B: multivariable Cox forest plot
# ============================================================

forest_df <- cox_source %>%
  mutate(
    term_label = dplyr::case_when(
      term == "four_gene_score_z" ~ "Four-gene score, per SD",
      term == "age" ~ "Age",
      grepl("^sex", term) ~ paste0("Sex: ", gsub("^sex", "", term)),
      grepl("^AJCC_stage", term) ~ paste0("AJCC stage: ", gsub("^AJCC_stage", "", term)),
      TRUE ~ term
    ),
    term_label = factor(term_label, levels = rev(term_label))
  )

forest_plot <- ggplot(
  forest_df,
  aes(x = HR, y = term_label)
) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_errorbarh(
    aes(xmin = lower_95_CI, xmax = upper_95_CI),
    height = 0.2,
    linewidth = 0.5
  ) +
  geom_point(size = 2) +
  scale_x_log10() +
  labs(
    x = "Hazard ratio, log scale",
    y = NULL
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(outdir, "Fig6B_TCGA_four_gene_score_multivariable_Cox.pdf"),
  plot = forest_plot,
  width = 7,
  height = 4.8,
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Fig6B_TCGA_four_gene_score_multivariable_Cox.png"),
  plot = forest_plot,
  width = 7,
  height = 4.8,
  dpi = 300,
  bg = "white"
)

# ============================================================
# Save session info
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  con = file.path(outdir, "07_LASSO_risk_score_sessionInfo.txt")
)

message("=== Script 07 complete ===")
message("Outputs written to: ", outdir)
message("Primary source data: 07_TCGA_four_gene_score_source_data.csv")
message("QC: 07_LASSO_risk_score_QC_check.csv")
