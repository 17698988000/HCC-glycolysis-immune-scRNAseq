# ============================================================
# 10_GSE14520_validation.R
#
# Purpose:
#   Independent validation of the TCGA-derived four-gene
#   glycolysis score in the GSE14520 HCC microarray cohort.
#
# Final manuscript alignment:
#   Dataset: GSE14520
#   Platform: GPL3921
#   Model genes: TPI1, ENO1, LDHA, SLC2A1
#   Score:
#     0.3041908*TPI1 +
#     0.9639654*ENO1 +
#     1.3404374*LDHA +
#     0.2424239*SLC2A1
#
# Expected manuscript-level results:
#   Tumor samples with complete survival information: n = 221
#   Multivariable Cox evaluable samples: n = 217
#   GSE14520 validation HR per SD: approximately 1.32
#
# Output policy:
#   1. Write source CSV.
#   2. Write QC CSV.
#   3. Generate PDF/PNG only if mandatory QC checks pass.
#
# ============================================================

options(stringsAsFactors = FALSE)

# -----------------------------
# 0. Package checks
# -----------------------------

required_pkgs <- c(
  "GEOquery",
  "Biobase",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "readr",
  "ggplot2",
  "survival",
  "survminer",
  "timeROC"
)

missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(missing_pkgs) > 0) {
  stop(
    "Missing required packages: ",
    paste(missing_pkgs, collapse = ", "),
    "\nInstall missing CRAN/Bioconductor packages before running this script."
  )
}

library(GEOquery)
library(Biobase)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(ggplot2)
library(survival)
library(survminer)
library(timeROC)

# Avoid common namespace masking
select      <- dplyr::select
mutate      <- dplyr::mutate
filter      <- dplyr::filter
arrange     <- dplyr::arrange
summarise   <- dplyr::summarise
group_by    <- dplyr::group_by
ungroup     <- dplyr::ungroup
left_join   <- dplyr::left_join
rename      <- dplyr::rename
bind_rows   <- dplyr::bind_rows
all_of      <- dplyr::all_of

set.seed(20260519)

# -----------------------------
# 1. User parameters
# -----------------------------

gse_id <- "GSE14520"
target_platform <- "GPL3921"

data_dir <- file.path("data", "GSE14520")
outdir <- file.path("results", "GSE14520_validation")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

supplement_file <- file.path(data_dir, "GSE14520_Extra_Supplement.txt.gz")
supplement_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz"

model_coef <- c(
  TPI1   = 0.3041908,
  ENO1   = 0.9639654,
  LDHA   = 1.3404374,
  SLC2A1 = 0.2424239
)

model_genes <- names(model_coef)

# Main median split rule.
# Consistent with the project convention: High = score > median; Low = remaining.
high_rule <- "greater_than_median"

# Manuscript targets used for QC reference.
expected_survival_n <- 221
expected_multivar_n <- 217
expected_hr_per_sd <- 1.32
expected_logrank_p <- 0.033
expected_multivar_p <- 0.0147
expected_cindex_clinical <- 0.652
expected_cindex_plus_score <- 0.672
expected_auc_1y <- 0.559
expected_auc_3y <- 0.605
expected_auc_5y <- 0.652

# -----------------------------
# 2. Helper functions
# -----------------------------

clean_key <- function(x) {
  x %>%
    tolower() %>%
    stringr::str_replace_all("[^a-z0-9]", "")
}

find_col <- function(df, patterns, label, required = TRUE) {
  cn <- names(df)
  key <- clean_key(cn)

  for (pat in patterns) {
    idx <- which(stringr::str_detect(key, pat))
    if (length(idx) > 0) {
      return(cn[idx[1]])
    }
  }

  if (required) {
    stop(
      "Could not identify required clinical column for: ", label,
      "\nAvailable columns:\n",
      paste(names(df), collapse = ", ")
    )
  }

  return(NA_character_)
}

parse_status_01 <- function(x) {
  s <- trimws(as.character(x))
  s_low <- tolower(s)

  out <- rep(NA_integer_, length(s))

  out[stringr::str_detect(s_low, "dead|deceased|death|died|event|yes|true")] <- 1L
  out[stringr::str_detect(s_low, "alive|living|censor|censored|no|false")] <- 0L

  numeric_x <- suppressWarnings(as.numeric(s))
  numeric_idx <- is.na(out) & !is.na(numeric_x)

  if (any(numeric_idx)) {
    vals <- sort(unique(numeric_x[numeric_idx]))

    if (all(vals %in% c(0, 1))) {
      out[numeric_idx] <- as.integer(numeric_x[numeric_idx])
    } else if (all(vals %in% c(1, 2))) {
      # Conservative fallback for 1/2 coding: use the larger value as event.
      out[numeric_idx] <- ifelse(numeric_x[numeric_idx] == max(vals), 1L, 0L)
    } else {
      warning(
        "Unrecognized numeric survival status coding: ",
        paste(vals, collapse = ", "),
        ". These values were set to NA."
      )
    }
  }

  out
}

parse_afp_group <- function(x) {
  s <- trimws(as.character(x))
  s_low <- tolower(s)
  num <- suppressWarnings(readr::parse_number(s))

  out <- rep(NA_character_, length(s))

  out[stringr::str_detect(s_low, ">|high|positive|elevated")] <- ">300"
  out[stringr::str_detect(s_low, "<=|≤|<|low|normal|negative")] <- "<=300"

  idx <- is.na(out) & !is.na(num)
  out[idx] <- ifelse(num[idx] > 300, ">300", "<=300")

  factor(out, levels = c("<=300", ">300"))
}

parse_yes_no <- function(x) {
  s <- trimws(as.character(x))
  s_low <- tolower(s)
  num <- suppressWarnings(as.numeric(s))

  out <- rep(NA_character_, length(s))

  out[stringr::str_detect(s_low, "yes|present|positive|true|multi|multiple")] <- "Yes"
  out[stringr::str_detect(s_low, "no|absent|negative|false|single|none")] <- "No"

  idx <- is.na(out) & !is.na(num)
  out[idx] <- ifelse(num[idx] > 0, "Yes", "No")

  factor(out, levels = c("No", "Yes"))
}

parse_size_group <- function(x) {
  s <- trimws(as.character(x))
  s_low <- tolower(s)
  num <- suppressWarnings(readr::parse_number(s))

  out <- rep(NA_character_, length(s))

  out[stringr::str_detect(s_low, "large|big|>|high")] <- "Large"
  out[stringr::str_detect(s_low, "small|<=|≤|<|low")] <- "Small"

  # Common HCC convention for large-vs-small tumor size is >5 cm.
  idx <- is.na(out) & !is.na(num)
  out[idx] <- ifelse(num[idx] > 5, "Large", "Small")

  factor(out, levels = c("Small", "Large"))
}

parse_multinodular_group <- function(x) {
  s <- trimws(as.character(x))
  s_low <- tolower(s)
  num <- suppressWarnings(readr::parse_number(s))

  out <- rep(NA_character_, length(s))

  out[stringr::str_detect(s_low, "multi|multiple|yes|present|positive|true")] <- "Yes"
  out[stringr::str_detect(s_low, "single|solitary|no|absent|negative|false")] <- "No"

  idx <- is.na(out) & !is.na(num)
  out[idx] <- ifelse(num[idx] > 1, "Yes", "No")

  factor(out, levels = c("No", "Yes"))
}

parse_tnm_group <- function(x) {
  s <- trimws(as.character(x))
  s_up <- toupper(s)
  num <- suppressWarnings(readr::parse_number(s))

  out <- rep(NA_character_, length(s))

  out[stringr::str_detect(s_up, "III|IV|3|4")] <- "TNM III-IV"
  out[stringr::str_detect(s_up, "I|II|1|2")] <- "TNM I-II"

  idx <- is.na(out) & !is.na(num)
  out[idx] <- ifelse(num[idx] >= 3, "TNM III-IV", "TNM I-II")

  factor(out, levels = c("TNM I-II", "TNM III-IV"))
}

extract_gene_symbols <- function(fdat, model_genes) {
  possible_symbol_cols <- c(
    "Gene Symbol",
    "Gene symbol",
    "GENE_SYMBOL",
    "SYMBOL",
    "Symbol",
    "gene_symbol",
    "Gene.symbol",
    "gene_assignment",
    "Gene Assignment"
  )

  symbol_col <- intersect(possible_symbol_cols, names(fdat))[1]

  if (is.na(symbol_col)) {
    candidate_cols <- names(fdat)[stringr::str_detect(clean_key(names(fdat)), "genesymbol|symbol|geneassignment")]
    if (length(candidate_cols) > 0) {
      symbol_col <- candidate_cols[1]
    }
  }

  if (is.na(symbol_col)) {
    stop(
      "Could not identify a gene-symbol annotation column in GPL annotation.\n",
      "Available fData columns:\n",
      paste(names(fdat), collapse = ", ")
    )
  }

  raw_symbol <- as.character(fdat[[symbol_col]])
  probe_id <- rownames(fdat)

  # Primary exact-style parsing for standard AnnotGPL "Gene Symbol" columns.
  parsed <- tibble(
    probe = probe_id,
    symbol_raw = raw_symbol
  ) %>%
    filter(!is.na(symbol_raw), symbol_raw != "", symbol_raw != "---") %>%
    mutate(
      symbol_raw = stringr::str_replace_all(symbol_raw, "\\s*///\\s*", ";"),
      symbol_raw = stringr::str_replace_all(symbol_raw, "\\s*//\\s*", ";")
    ) %>%
    tidyr::separate_rows(symbol_raw, sep = "\\s*;\\s*|\\s*,\\s*") %>%
    mutate(
      symbol = stringr::str_trim(symbol_raw),
      symbol = stringr::str_replace(symbol, "\\s+.*$", "")
    ) %>%
    filter(symbol %in% model_genes) %>%
    distinct(probe, symbol)

  # Fallback for nonstandard gene-assignment strings.
  if (nrow(parsed) == 0) {
    parsed <- bind_rows(lapply(model_genes, function(g) {
      idx <- stringr::str_detect(
        raw_symbol,
        paste0("(^|[^A-Za-z0-9_.-])", g, "([^A-Za-z0-9_.-]|$)")
      )
      tibble(
        probe = probe_id[idx],
        symbol = g
      )
    })) %>%
      distinct(probe, symbol)
  }

  attr(parsed, "symbol_col") <- symbol_col
  parsed
}

format_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, formatC(p, format = "e", digits = 2), signif(p, 3))
  )
}

qc_row <- function(check, observed, expected, pass, severity = "ERROR", notes = "") {
  tibble(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    pass = as.logical(pass),
    severity = severity,
    notes = notes
  )
}

get_cindex <- function(fit) {
  as.numeric(summary(fit)$concordance[1])
}

# -----------------------------
# 3. Download/read GSE14520
# -----------------------------

cat("Loading GSE14520 expression data from GEO...\n")

gse <- GEOquery::getGEO(
  gse_id,
  GSEMatrix = TRUE,
  AnnotGPL = TRUE,
  getGPL = TRUE,
  destdir = data_dir
)

if (!is.list(gse)) {
  gse <- list(gse)
}

platforms <- vapply(gse, Biobase::annotation, character(1))
cat("Available platforms:\n")
print(platforms)

platform_idx <- which(platforms == target_platform)

if (length(platform_idx) == 0) {
  stop(
    "Target platform ", target_platform, " was not found in downloaded GSE object.\n",
    "Available platforms: ", paste(platforms, collapse = ", ")
  )
}

eset <- gse[[platform_idx[1]]]
expr_raw <- Biobase::exprs(eset)
fdat <- Biobase::fData(eset)
pdat <- Biobase::pData(eset) %>%
  tibble::rownames_to_column("sample")

cat("Selected platform:", Biobase::annotation(eset), "\n")
cat("Expression matrix dimensions:", paste(dim(expr_raw), collapse = " x "), "\n")

# GEO series matrix values are expected to be processed microarray expression.
# If a user accidentally supplies an unlogged matrix, transform defensively.
expr_q99 <- as.numeric(stats::quantile(expr_raw, probs = 0.99, na.rm = TRUE))
expr_transformed <- FALSE

expr_use <- expr_raw
if (expr_q99 > 100) {
  expr_use <- log2(expr_raw + 1)
  expr_transformed <- TRUE
  warning("Expression values appeared unlogged; applied log2(x + 1).")
}

# -----------------------------
# 4. Download/read clinical supplement
# -----------------------------

if (!file.exists(supplement_file)) {
  cat("Downloading GSE14520 clinical supplement...\n")
  download.file(
    supplement_url,
    destfile = supplement_file,
    mode = "wb",
    quiet = FALSE
  )
}

if (!file.exists(supplement_file)) {
  stop("Clinical supplement file was not found and could not be downloaded: ", supplement_file)
}

clin_raw <- readr::read_tsv(
  supplement_file,
  show_col_types = FALSE,
  progress = FALSE
)

cat("Clinical supplement dimensions:", paste(dim(clin_raw), collapse = " x "), "\n")
cat("Clinical supplement columns:\n")
print(names(clin_raw))

sample_col <- find_col(
  clin_raw,
  patterns = c("^affygsm$", "affygsm", "^gsm$", "geoaccession", "sample"),
  label = "sample/GSM identifier"
)

tissue_col <- find_col(
  clin_raw,
  patterns = c("^tissuetype$", "tissuetype", "^tissue$", "sampletype"),
  label = "tissue type",
  required = FALSE
)

os_time_col <- find_col(
  clin_raw,
  patterns = c("survivalmonths", "survivalmonth", "osmonths", "overallsurvivalmonths", "ostime"),
  label = "overall survival time in months"
)

os_status_col <- find_col(
  clin_raw,
  patterns = c("survivalstatus", "osstatus", "overallstatus", "vitalstatus", "death"),
  label = "overall survival status"
)

afp_col <- find_col(
  clin_raw,
  patterns = c("^afp$", "afpngml", "afplevel", "alphafetoprotein"),
  label = "AFP",
  required = FALSE
)

cirrhosis_col <- find_col(
  clin_raw,
  patterns = c("cirrhosis", "cirrhotic"),
  label = "cirrhosis",
  required = FALSE
)

tumor_size_col <- find_col(
  clin_raw,
  patterns = c("maintumorsize", "tumorsize", "sizeofmaintumor", "maximaldiameter", "diameter"),
  label = "main tumor size",
  required = FALSE
)

multinodular_col <- find_col(
  clin_raw,
  patterns = c("multinodular", "multinodule", "tumornumber", "nodular", "nodule"),
  label = "multinodular disease",
  required = FALSE
)

tnm_col <- find_col(
  clin_raw,
  patterns = c("^tnm$", "tnmstage", "stage"),
  label = "TNM stage",
  required = FALSE
)

clinical_df <- clin_raw %>%
  transmute(
    sample = as.character(.data[[sample_col]]),
    tissue_type_clinical = if (!is.na(tissue_col)) as.character(.data[[tissue_col]]) else NA_character_,
    OS_time_months = suppressWarnings(as.numeric(.data[[os_time_col]])),
    OS_status = parse_status_01(.data[[os_status_col]]),
    AFP_group = if (!is.na(afp_col)) parse_afp_group(.data[[afp_col]]) else factor(NA_character_, levels = c("<=300", ">300")),
    cirrhosis_group = if (!is.na(cirrhosis_col)) parse_yes_no(.data[[cirrhosis_col]]) else factor(NA_character_, levels = c("No", "Yes")),
    main_tumor_size_group = if (!is.na(tumor_size_col)) parse_size_group(.data[[tumor_size_col]]) else factor(NA_character_, levels = c("Small", "Large")),
    multinodular_group = if (!is.na(multinodular_col)) parse_multinodular_group(.data[[multinodular_col]]) else factor(NA_character_, levels = c("No", "Yes")),
    TNM_group = if (!is.na(tnm_col)) parse_tnm_group(.data[[tnm_col]]) else factor(NA_character_, levels = c("TNM I-II", "TNM III-IV"))
  ) %>%
  distinct(sample, .keep_all = TRUE)

geo_sample_df <- pdat %>%
  mutate(
    geo_title = if ("title" %in% names(.)) as.character(.data$title) else NA_character_,
    tissue_type_geo = dplyr::case_when(
      stringr::str_detect(tolower(geo_title), "non[- ]?tumor|normal") ~ "Non-Tumor",
      stringr::str_detect(tolower(geo_title), "tumor") ~ "Tumor",
      TRUE ~ NA_character_
    )
  ) %>%
  select(sample, geo_title, tissue_type_geo)

clinical_df <- clinical_df %>%
  left_join(geo_sample_df, by = "sample") %>%
  mutate(
    tissue_type = dplyr::coalesce(tissue_type_clinical, tissue_type_geo),
    is_tumor = stringr::str_detect(tolower(tissue_type), "tumor") &
      !stringr::str_detect(tolower(tissue_type), "non")
  )

# -----------------------------
# 5. Probe-to-gene collapse
# -----------------------------

cat("Collapsing probes to model genes...\n")

probe_symbol <- extract_gene_symbols(fdat, model_genes)
symbol_col_used <- attr(probe_symbol, "symbol_col")

cat("Annotation symbol column used:", symbol_col_used, "\n")
cat("Probe coverage by model gene:\n")
print(table(probe_symbol$symbol))

missing_model_genes <- setdiff(model_genes, unique(probe_symbol$symbol))

if (length(missing_model_genes) > 0) {
  warning("Missing model genes in platform annotation: ", paste(missing_model_genes, collapse = ", "))
}

expr_tbl <- as.data.frame(expr_use, check.names = FALSE) %>%
  tibble::rownames_to_column("probe")

gene_expr_long <- expr_tbl %>%
  inner_join(probe_symbol, by = "probe") %>%
  tidyr::pivot_longer(
    cols = all_of(colnames(expr_use)),
    names_to = "sample",
    values_to = "expression"
  ) %>%
  group_by(symbol, sample) %>%
  summarise(
    expression = mean(expression, na.rm = TRUE),
    n_probes = dplyr::n(),
    .groups = "drop"
  )

gene_expr_wide <- gene_expr_long %>%
  select(symbol, sample, expression) %>%
  tidyr::pivot_wider(
    names_from = symbol,
    values_from = expression
  )

probe_count_wide <- gene_expr_long %>%
  select(symbol, sample, n_probes) %>%
  tidyr::pivot_wider(
    names_from = symbol,
    values_from = n_probes,
    names_prefix = "n_probe_"
  )

for (g in model_genes) {
  if (!g %in% names(gene_expr_wide)) {
    gene_expr_wide[[g]] <- NA_real_
  }
  probe_col <- paste0("n_probe_", g)
  if (!probe_col %in% names(probe_count_wide)) {
    probe_count_wide[[probe_col]] <- NA_integer_
  }
}

gene_expr_wide <- gene_expr_wide %>%
  select(sample, all_of(model_genes))

probe_count_wide <- probe_count_wide %>%
  select(sample, all_of(paste0("n_probe_", model_genes)))

# -----------------------------
# 6. Build analysis table
# -----------------------------

analysis_df <- clinical_df %>%
  filter(is_tumor) %>%
  inner_join(gene_expr_wide, by = "sample") %>%
  left_join(probe_count_wide, by = "sample") %>%
  mutate(
    four_gene_score =
      model_coef["TPI1"]   * TPI1 +
      model_coef["ENO1"]   * ENO1 +
      model_coef["LDHA"]   * LDHA +
      model_coef["SLC2A1"] * SLC2A1
  ) %>%
  filter(
    !is.na(OS_time_months),
    !is.na(OS_status),
    !is.na(four_gene_score)
  ) %>%
  mutate(
    score_z = as.numeric(scale(four_gene_score)),
    score_median = median(four_gene_score, na.rm = TRUE),
    score_group = ifelse(four_gene_score > score_median, "High", "Low"),
    score_group = factor(score_group, levels = c("Low", "High"))
  ) %>%
  arrange(sample)

source_file <- file.path(outdir, "FigS17_GSE14520_four_gene_score_source_data.csv")
readr::write_csv(analysis_df, source_file)

cat("Source data written to:", source_file, "\n")
cat("Survival-complete tumor sample count:", nrow(analysis_df), "\n")
cat("Score group counts:\n")
print(table(analysis_df$score_group))

# -----------------------------
# 7. Survival analyses
# -----------------------------

km_fit <- survival::survfit(
  survival::Surv(OS_time_months, OS_status) ~ score_group,
  data = analysis_df
)

logrank <- survival::survdiff(
  survival::Surv(OS_time_months, OS_status) ~ score_group,
  data = analysis_df
)

logrank_p <- 1 - stats::pchisq(logrank$chisq, df = 1)

cox_group <- survival::coxph(
  survival::Surv(OS_time_months, OS_status) ~ score_group,
  data = analysis_df
)

cox_score_uni <- survival::coxph(
  survival::Surv(OS_time_months, OS_status) ~ score_z,
  data = analysis_df
)

cox_group_summary <- summary(cox_group)
cox_score_uni_summary <- summary(cox_score_uni)

# Multivariable model specified in final manuscript:
# AFP level, cirrhosis status, main tumor size, multinodular disease.
mv_df <- analysis_df %>%
  filter(
    !is.na(AFP_group),
    !is.na(cirrhosis_group),
    !is.na(main_tumor_size_group),
    !is.na(multinodular_group),
    !is.na(score_z)
  ) %>%
  mutate(
    AFP_group = stats::relevel(AFP_group, ref = "<=300"),
    cirrhosis_group = stats::relevel(cirrhosis_group, ref = "No"),
    main_tumor_size_group = stats::relevel(main_tumor_size_group, ref = "Small"),
    multinodular_group = stats::relevel(multinodular_group, ref = "No")
  )

clinical_formula <- survival::Surv(OS_time_months, OS_status) ~
  AFP_group +
  cirrhosis_group +
  main_tumor_size_group +
  multinodular_group

plus_score_formula <- survival::Surv(OS_time_months, OS_status) ~
  AFP_group +
  cirrhosis_group +
  main_tumor_size_group +
  multinodular_group +
  score_z

fit_clinical <- survival::coxph(
  clinical_formula,
  data = mv_df,
  x = TRUE,
  y = TRUE
)

fit_plus_score <- survival::coxph(
  plus_score_formula,
  data = mv_df,
  x = TRUE,
  y = TRUE
)

fit_plus_summary <- summary(fit_plus_score)

score_row <- which(rownames(fit_plus_summary$coefficients) == "score_z")

if (length(score_row) != 1) {
  stop("Could not locate score_z in multivariable Cox summary.")
}

hr_per_sd <- as.numeric(fit_plus_summary$conf.int[score_row, "exp(coef)"])
hr_lower <- as.numeric(fit_plus_summary$conf.int[score_row, "lower .95"])
hr_upper <- as.numeric(fit_plus_summary$conf.int[score_row, "upper .95"])
hr_p <- as.numeric(fit_plus_summary$coefficients[score_row, "Pr(>|z|)"])

aic_clinical <- stats::AIC(fit_clinical)
aic_plus_score <- stats::AIC(fit_plus_score)

cindex_clinical <- get_cindex(fit_clinical)
cindex_plus_score <- get_cindex(fit_plus_score)

lrt_tab <- stats::anova(fit_clinical, fit_plus_score, test = "LRT")
lrt_chisq <- as.numeric(lrt_tab$Chisq[2])
lrt_p <- as.numeric(lrt_tab$`Pr(>|Chi|)`[2])

ph_test <- survival::cox.zph(fit_plus_score)

# -----------------------------
# 8. Time-dependent ROC
# -----------------------------

roc_times <- c(12, 36, 60)

time_roc <- timeROC::timeROC(
  T = analysis_df$OS_time_months,
  delta = analysis_df$OS_status,
  marker = analysis_df$four_gene_score,
  cause = 1,
  times = roc_times,
  iid = TRUE
)

time_roc_df <- tibble(
  time_months = roc_times,
  auc = as.numeric(time_roc$AUC)
)

time_roc_file <- file.path(outdir, "TableS4_GSE14520_time_dependent_AUC.csv")
readr::write_csv(time_roc_df, time_roc_file)

# -----------------------------
# 9. TNM-stratified exploratory analyses
# -----------------------------

tnm_df <- analysis_df %>%
  filter(!is.na(TNM_group)) %>%
  mutate(TNM_group = factor(TNM_group, levels = c("TNM I-II", "TNM III-IV")))

tnm_results <- tibble()

if (nrow(tnm_df) > 0 && length(unique(stats::na.omit(tnm_df$TNM_group))) >= 1) {
  tnm_results <- bind_rows(lapply(levels(tnm_df$TNM_group), function(grp) {
    dat <- tnm_df %>% filter(TNM_group == grp)

    if (nrow(dat) < 10 ||
        length(unique(dat$score_group)) < 2 ||
        sum(dat$OS_status == 1, na.rm = TRUE) < 2) {
      return(tibble(
        TNM_group = grp,
        n = nrow(dat),
        events = sum(dat$OS_status == 1, na.rm = TRUE),
        HR_high_vs_low = NA_real_,
        HR_lower = NA_real_,
        HR_upper = NA_real_,
        logrank_p = NA_real_,
        notes = "Insufficient events or group variation"
      ))
    }

    lr <- survival::survdiff(
      survival::Surv(OS_time_months, OS_status) ~ score_group,
      data = dat
    )
    lr_p <- 1 - stats::pchisq(lr$chisq, df = 1)

    cx <- survival::coxph(
      survival::Surv(OS_time_months, OS_status) ~ score_group,
      data = dat
    )
    cx_sum <- summary(cx)

    tibble(
      TNM_group = grp,
      n = nrow(dat),
      events = sum(dat$OS_status == 1, na.rm = TRUE),
      HR_high_vs_low = as.numeric(cx_sum$conf.int[1, "exp(coef)"]),
      HR_lower = as.numeric(cx_sum$conf.int[1, "lower .95"]),
      HR_upper = as.numeric(cx_sum$conf.int[1, "upper .95"]),
      logrank_p = lr_p,
      notes = ""
    )
  }))
}

tnm_file <- file.path(outdir, "TableS4_GSE14520_TNM_stratified_survival.csv")
readr::write_csv(tnm_results, tnm_file)

# -----------------------------
# 10. Model metrics table
# -----------------------------

model_metrics <- tibble(
  metric = c(
    "survival_complete_tumor_n",
    "multivariable_complete_n",
    "median_score",
    "score_group_low_n",
    "score_group_high_n",
    "logrank_p",
    "univariate_score_HR_per_SD",
    "univariate_score_p",
    "multivariable_score_HR_per_SD",
    "multivariable_score_HR_lower_95",
    "multivariable_score_HR_upper_95",
    "multivariable_score_p",
    "AIC_clinical_only",
    "AIC_clinical_plus_score",
    "C_index_clinical_only",
    "C_index_clinical_plus_score",
    "LRT_chisq",
    "LRT_p",
    "timeROC_AUC_12_months",
    "timeROC_AUC_36_months",
    "timeROC_AUC_60_months"
  ),
  value = c(
    nrow(analysis_df),
    nrow(mv_df),
    unique(analysis_df$score_median)[1],
    as.integer(table(analysis_df$score_group)["Low"]),
    as.integer(table(analysis_df$score_group)["High"]),
    logrank_p,
    as.numeric(cox_score_uni_summary$conf.int[1, "exp(coef)"]),
    as.numeric(cox_score_uni_summary$coefficients[1, "Pr(>|z|)"]),
    hr_per_sd,
    hr_lower,
    hr_upper,
    hr_p,
    aic_clinical,
    aic_plus_score,
    cindex_clinical,
    cindex_plus_score,
    lrt_chisq,
    lrt_p,
    time_roc_df$auc[time_roc_df$time_months == 12],
    time_roc_df$auc[time_roc_df$time_months == 36],
    time_roc_df$auc[time_roc_df$time_months == 60]
  )
)

metrics_file <- file.path(outdir, "TableS4_GSE14520_model_metrics.csv")
readr::write_csv(model_metrics, metrics_file)

ph_file <- file.path(outdir, "TableS4_GSE14520_PH_test.csv")
ph_df <- as.data.frame(ph_test$table) %>%
  tibble::rownames_to_column("term")
readr::write_csv(ph_df, ph_file)

# -----------------------------
# 11. QC table
# -----------------------------

gene_presence <- setNames(model_genes %in% unique(probe_symbol$symbol), model_genes)
all_genes_present <- all(gene_presence)

qc <- bind_rows(
  qc_row(
    "selected_platform",
    Biobase::annotation(eset),
    target_platform,
    Biobase::annotation(eset) == target_platform
  ),
  qc_row(
    "clinical_supplement_file_exists",
    file.exists(supplement_file),
    TRUE,
    file.exists(supplement_file)
  ),
  qc_row(
    "expression_transformed_log2_if_needed",
    expr_transformed,
    "FALSE expected for processed GEO series matrix; TRUE allowed if raw-like scale detected",
    TRUE,
    severity = "INFO"
  ),
  qc_row(
    "model_gene_presence",
    paste(names(gene_presence), gene_presence, sep = "=", collapse = "; "),
    "All four genes present: TPI1, ENO1, LDHA, SLC2A1",
    all_genes_present
  ),
  qc_row(
    "survival_complete_tumor_n",
    nrow(analysis_df),
    expected_survival_n,
    nrow(analysis_df) == expected_survival_n
  ),
  qc_row(
    "multivariable_complete_n",
    nrow(mv_df),
    expected_multivar_n,
    nrow(mv_df) == expected_multivar_n
  ),
  qc_row(
    "missing_four_gene_score_n",
    sum(is.na(analysis_df$four_gene_score)),
    0,
    sum(is.na(analysis_df$four_gene_score)) == 0
  ),
  qc_row(
    "score_group_counts",
    paste(names(table(analysis_df$score_group)), as.integer(table(analysis_df$score_group)), sep = "=", collapse = "; "),
    "Median split; High = score > median; Low = remaining",
    length(unique(analysis_df$score_group)) == 2
  ),
  qc_row(
    "logrank_p_directional",
    signif(logrank_p, 4),
    paste0("Manuscript target approximately ", expected_logrank_p, "; must be < 0.05"),
    !is.na(logrank_p) && logrank_p < 0.05
  ),
  qc_row(
    "multivariable_HR_per_SD_directional",
    signif(hr_per_sd, 4),
    paste0("Manuscript target approximately ", expected_hr_per_sd, "; must be > 1"),
    !is.na(hr_per_sd) && hr_per_sd > 1
  ),
  qc_row(
    "multivariable_score_p_directional",
    signif(hr_p, 4),
    paste0("Manuscript target approximately ", expected_multivar_p, "; must be < 0.05"),
    !is.na(hr_p) && hr_p < 0.05
  ),
  qc_row(
    "clinical_plus_score_AIC_improvement",
    paste0("clinical=", round(aic_clinical, 3), "; plus_score=", round(aic_plus_score, 3)),
    "AIC should decrease after adding score",
    aic_plus_score < aic_clinical
  ),
  qc_row(
    "clinical_plus_score_C_index_improvement",
    paste0("clinical=", round(cindex_clinical, 3), "; plus_score=", round(cindex_plus_score, 3)),
    "C-index should increase after adding score",
    cindex_plus_score > cindex_clinical
  ),
  qc_row(
    "LRT_p",
    signif(lrt_p, 4),
    "Likelihood-ratio p should be < 0.05",
    !is.na(lrt_p) && lrt_p < 0.05
  ),
  qc_row(
    "reference_logrank_p",
    signif(logrank_p, 4),
    expected_logrank_p,
    abs(logrank_p - expected_logrank_p) <= 0.02,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  ),
  qc_row(
    "reference_multivariable_HR_per_SD",
    signif(hr_per_sd, 4),
    expected_hr_per_sd,
    abs(hr_per_sd - expected_hr_per_sd) <= 0.15,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  ),
  qc_row(
    "reference_C_index_clinical_only",
    signif(cindex_clinical, 4),
    expected_cindex_clinical,
    abs(cindex_clinical - expected_cindex_clinical) <= 0.05,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  ),
  qc_row(
    "reference_C_index_clinical_plus_score",
    signif(cindex_plus_score, 4),
    expected_cindex_plus_score,
    abs(cindex_plus_score - expected_cindex_plus_score) <= 0.05,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  ),
  qc_row(
    "reference_timeROC_AUC_12_months",
    signif(time_roc_df$auc[time_roc_df$time_months == 12], 4),
    expected_auc_1y,
    abs(time_roc_df$auc[time_roc_df$time_months == 12] - expected_auc_1y) <= 0.08,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  ),
  qc_row(
    "reference_timeROC_AUC_36_months",
    signif(time_roc_df$auc[time_roc_df$time_months == 36], 4),
    expected_auc_3y,
    abs(time_roc_df$auc[time_roc_df$time_months == 36] - expected_auc_3y) <= 0.08,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  ),
  qc_row(
    "reference_timeROC_AUC_60_months",
    signif(time_roc_df$auc[time_roc_df$time_months == 60], 4),
    expected_auc_5y,
    abs(time_roc_df$auc[time_roc_df$time_months == 60] - expected_auc_5y) <= 0.08,
    severity = "CHECK",
    notes = "Reference tolerance check; not used to block figure output."
  )
)

qc_file <- file.path(outdir, "FigS17_GSE14520_QC_check.csv")
readr::write_csv(qc, qc_file)

cat("QC table written to:", qc_file, "\n")
print(qc)

mandatory_pass <- all(qc$pass[qc$severity == "ERROR"])

if (!mandatory_pass) {
  stop(
    "Mandatory QC failed. Source and QC CSV files were written, but PDF/PNG were not generated. ",
    "Inspect: ", qc_file
  )
}

cat("Mandatory QC PASS. Proceeding to PDF/PNG generation.\n")

# -----------------------------
# 12. Plot Supplementary Figure S17
# -----------------------------

km_label <- paste0(
  "Log-rank p = ", format_p(logrank_p),
  "\nHR per SD = ", signif(hr_per_sd, 3),
  " (95% CI ", signif(hr_lower, 3), "-", signif(hr_upper, 3), ")"
)

p_km <- survminer::ggsurvplot(
  km_fit,
  data = analysis_df,
  pval = FALSE,
  risk.table = TRUE,
  risk.table.height = 0.25,
  conf.int = FALSE,
  censor = TRUE,
  palette = c("#4DBBD5", "#E64B35"),
  legend.title = "",
  legend.labs = c("Low four-gene score", "High four-gene score"),
  xlab = "Time (months)",
  ylab = "Overall survival probability",
  break.time.by = 20,
  ggtheme = ggplot2::theme_classic(base_size = 12),
  risk.table.y.text = FALSE,
  risk.table.title = "Number at risk"
)

p_km$plot <- p_km$plot +
  ggplot2::annotate(
    "text",
    x = 0,
    y = 0.12,
    label = km_label,
    hjust = 0,
    vjust = 0,
    size = 3.6
  ) +
  ggplot2::theme(
    legend.position = c(0.78, 0.86),
    legend.background = ggplot2::element_blank(),
    plot.title = ggplot2::element_blank()
  )

fig_pdf <- file.path(outdir, "FigS17_GSE14520_KM_four_gene_score.pdf")
fig_png <- file.path(outdir, "FigS17_GSE14520_KM_four_gene_score.png")

ggplot2::ggsave(
  filename = fig_pdf,
  plot = print(p_km),
  width = 7.2,
  height = 7.2,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = fig_png,
  plot = print(p_km),
  width = 7.2,
  height = 7.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

# -----------------------------
# 13. Optional ROC plot
# -----------------------------

roc_plot_df <- time_roc_df %>%
  mutate(
    time_label = dplyr::case_when(
      time_months == 12 ~ "1 year",
      time_months == 36 ~ "3 years",
      time_months == 60 ~ "5 years",
      TRUE ~ paste0(time_months, " months")
    )
  )

p_auc <- ggplot2::ggplot(
  roc_plot_df,
  ggplot2::aes(x = factor(time_label, levels = c("1 year", "3 years", "5 years")), y = auc)
) +
  ggplot2::geom_col(width = 0.65) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.3f", auc)),
    vjust = -0.4,
    size = 3.8
  ) +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::labs(
    x = NULL,
    y = "Time-dependent AUC"
  ) +
  ggplot2::theme_classic(base_size = 12)

auc_pdf <- file.path(outdir, "GSE14520_time_dependent_AUC_four_gene_score.pdf")
auc_png <- file.path(outdir, "GSE14520_time_dependent_AUC_four_gene_score.png")

ggplot2::ggsave(
  filename = auc_pdf,
  plot = p_auc,
  width = 4.8,
  height = 4.2,
  units = "in",
  bg = "white"
)

ggplot2::ggsave(
  filename = auc_png,
  plot = p_auc,
  width = 4.8,
  height = 4.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

# -----------------------------
# 14. Console summary
# -----------------------------

cat("\nGSE14520 four-gene validation complete.\n")
cat("Source CSV: ", source_file, "\n")
cat("QC CSV: ", qc_file, "\n")
cat("Metrics CSV: ", metrics_file, "\n")
cat("TimeROC CSV: ", time_roc_file, "\n")
cat("TNM CSV: ", tnm_file, "\n")
cat("PH test CSV: ", ph_file, "\n")
cat("KM PDF: ", fig_pdf, "\n")
cat("KM PNG: ", fig_png, "\n")
cat("AUC PDF: ", auc_pdf, "\n")
cat("AUC PNG: ", auc_png, "\n")

cat("\nKey metrics:\n")
cat("n survival-complete tumor samples: ", nrow(analysis_df), "\n")
cat("n multivariable complete samples: ", nrow(mv_df), "\n")
cat("log-rank p: ", logrank_p, "\n")
cat("HR per SD: ", hr_per_sd, " [", hr_lower, ", ", hr_upper, "]\n", sep = "")
cat("multivariable p: ", hr_p, "\n")
cat("AIC clinical only: ", aic_clinical, "\n")
cat("AIC clinical + score: ", aic_plus_score, "\n")
cat("C-index clinical only: ", cindex_clinical, "\n")
cat("C-index clinical + score: ", cindex_plus_score, "\n")
cat("LRT p: ", lrt_p, "\n")
cat("AUC 12/36/60 months: ", paste(round(time_roc_df$auc, 3), collapse = " / "), "\n")
