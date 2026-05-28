# ============================================================
# 12_drug_repurposing.R
#
# Current purpose:
#   Exploratory GSE235863 anti-PD-1 plus lenvatinib non-response
#   association analysis using the locked TCGA-derived four-gene
#   glycolysis score.
#
# This file replaces the legacy drug-target Fisher enrichment script.
# It does not perform drug-repurposing prioritization and does not
# make treatment-selection claims.
#
# Locked scientific settings:
#   Cohort: GSE235863, anti-PD-1 plus lenvatinib, HCC
#   Expected evaluable samples: n = 15
#   Responders: n = 11
#   Non-responders: n = 4
#   Expression scale used for scoring: log2(TPM + 1)
#   Positive ROC class: non-responder
#   Score:
#     0.3041908*TPI1 + 0.9639654*ENO1 +
#     1.3404374*LDHA + 0.2424239*SLC2A1
#   Expected exploratory result:
#     ROC AUC approximately 0.932
#     95% CI approximately 0.773-1.000
#     Median-based stratification with median sample included in
#     High group places all four non-responders in High group
#     Fisher exact p approximately 0.077
#
# Mandatory output order:
#   1. Source CSV files
#   2. QC CSV file
#   3. PDF/PNG figures only if QC PASS
#
# Interpretation:
#   Hypothesis-generating / exploratory only.
#   Not a formal predictive model.
#   Not a clinical treatment-selection assay.
# ============================================================

options(stringsAsFactors = FALSE)

# ----------------------------
# Package handling
# ----------------------------
required_packages <- c(
  "dplyr", "tidyr", "readr", "stringr", "tibble",
  "ggplot2", "pROC", "scales"
)

optional_bioc_packages <- c("GEOquery", "Biobase")

load_or_stop <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required but not installed.\n",
      "Install it before running this script.\n",
      "For CRAN packages use install.packages('", pkg, "').\n",
      "For Bioconductor packages use BiocManager::install('", pkg, "').",
      call. = FALSE
    )
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

invisible(lapply(required_packages, load_or_stop))

# Bioconductor packages are used for automatic GEO download / metadata retrieval.
# The script can still run from local files if GEOquery is unavailable.
geoquery_available <- requireNamespace("GEOquery", quietly = TRUE)
biobase_available <- requireNamespace("Biobase", quietly = TRUE)

# Avoid select() masking in interactive sessions.
select       <- dplyr::select
mutate       <- dplyr::mutate
filter       <- dplyr::filter
arrange      <- dplyr::arrange
summarise    <- dplyr::summarise
group_by     <- dplyr::group_by
ungroup      <- dplyr::ungroup
left_join    <- dplyr::left_join
rename       <- dplyr::rename
bind_rows    <- dplyr::bind_rows
n_distinct   <- dplyr::n_distinct
all_of       <- dplyr::all_of

# ----------------------------
# User-editable inputs
# ----------------------------
gse_id <- "GSE235863"

outdir <- "results_12_GSE235863_nonresponse"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# If these files exist, they are used preferentially.
# Expression file should contain TPM values, with genes as rows or in long format.
# Metadata file is optional and should contain at least:
#   sample_id,response
# where response can be R, responder, CR, PR, NR, non-responder, SD, or PD.
local_expression_file <- "GSE235863_TPM_matrix.csv"
local_metadata_file   <- "GSE235863_sample_metadata.csv"

# If expression values are already log2(TPM + 1), set this to TRUE.
# Locked analysis assumes log2(TPM + 1) is used for scoring.
input_expression_already_log2_tpm <- FALSE

# Automatic GEO download is attempted only if local_expression_file is absent.
allow_geo_download <- TRUE

# Expected locked QC values.
expected_n_total <- 15L
expected_n_responder <- 11L
expected_n_nonresponder <- 4L
expected_auc <- 0.932
expected_ci_low <- 0.773
expected_ci_high <- 1.000
expected_fisher_p <- 0.077

auc_tolerance <- 0.010
ci_tolerance <- 0.030
fisher_tolerance <- 0.020

# Locked four-gene score coefficients.
score_coefficients <- c(
  TPI1   = 0.3041908,
  ENO1   = 0.9639654,
  LDHA   = 1.3404374,
  SLC2A1 = 0.2424239
)
required_genes <- names(score_coefficients)

# ----------------------------
# Utility functions
# ----------------------------
clean_gene_symbol <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace(x, "\\.\\d+$", "")
  x <- stringr::str_replace_all(x, "\"", "")
  x <- stringr::str_trim(x)
  toupper(x)
}

clean_sample_id <- function(x) {
  x <- as.character(x)
  x <- stringr::str_replace_all(x, "\"", "")
  x <- stringr::str_trim(x)
  x
}

extract_patient_id <- function(x) {
  x <- as.character(x)
  hit <- stringr::str_match(x, stringr::regex("\\b(P\\d{1,3})\\b", ignore_case = TRUE))[, 2]
  toupper(hit)
}

standardize_response <- function(x, response_field_only = FALSE) {
  x0 <- as.character(x)
  x1 <- stringr::str_to_lower(x0)

  # Non-response is evaluated first to avoid matching "responder" inside
  # "non-responder".
  is_nonresponse <- stringr::str_detect(
    x1,
    "non[-_ ]?responder|\\bnr\\b|progressive disease|stable disease"
  )

  # In explicit response fields only, SD and PD are safe response labels.
  # In broader sample text, PD also appears in anti-PD-1 and is therefore
  # not used unless response_field_only = TRUE.
  if (response_field_only) {
    is_nonresponse <- is_nonresponse |
      stringr::str_detect(x1, "\\bpd\\b|\\bsd\\b")
  } else {
    is_nonresponse <- is_nonresponse |
      stringr::str_detect(x1, "\\bnr/")
  }

  is_response <- stringr::str_detect(
    x1,
    "\\bresponder\\b|\\br\\b|partial response|complete response|\\bpr\\b|\\bcr\\b"
  )
  is_response <- is_response & !is_nonresponse

  dplyr::case_when(
    is_nonresponse ~ "NonResponder",
    is_response ~ "Responder",
    TRUE ~ NA_character_
  )
}

extract_response_field <- function(text) {
  text <- as.character(text)
  # GEO characteristics often contain strings such as:
  # "response: PD, non-responder"
  m <- stringr::str_match(
    text,
    stringr::regex("response\\s*[:=]\\s*([^;|\\n]+)", ignore_case = TRUE)
  )[, 2]
  m
}

read_table_flexible <- function(path) {
  message("Reading candidate table: ", path)
  ext <- stringr::str_to_lower(basename(path))

  read_attempts <- list(
    function() readr::read_csv(path, show_col_types = FALSE, progress = FALSE),
    function() readr::read_tsv(path, show_col_types = FALSE, progress = FALSE),
    function() readr::read_delim(path, delim = "\t", show_col_types = FALSE, progress = FALSE),
    function() readr::read_delim(path, delim = ",", show_col_types = FALSE, progress = FALSE),
    function() read.delim(path, check.names = FALSE)
  )

  for (fn in read_attempts) {
    obj <- tryCatch(fn(), error = function(e) NULL)
    if (!is.null(obj) && nrow(obj) > 0 && ncol(obj) > 1) {
      obj <- as.data.frame(obj, check.names = FALSE)
      names(obj) <- make.unique(names(obj))
      return(obj)
    }
  }

  stop("Could not parse table: ", path, call. = FALSE)
}

numeric_column_score <- function(x) {
  y <- suppressWarnings(as.numeric(as.character(x)))
  sum(!is.na(y))
}

detect_expression_matrix <- function(df, source_file) {
  if (nrow(df) == 0 || ncol(df) < 3) return(NULL)

  cn <- names(df)
  cn_lower <- stringr::str_to_lower(cn)

  gene_col_candidates <- cn[vapply(df, function(col) {
    vals <- clean_gene_symbol(col)
    sum(vals %in% required_genes, na.rm = TRUE) >= 2
  }, logical(1))]

  # Long-format support: gene, sample, value/TPM columns.
  if (length(gene_col_candidates) > 0) {
    gene_col <- gene_col_candidates[1]

    sample_col_candidates <- cn[
      stringr::str_detect(cn_lower, "sample|gsm|patient|library|id") &
        cn != gene_col
    ]

    value_col_candidates <- cn[
      stringr::str_detect(cn_lower, "tpm|expression|expr|value|abundance") &
        cn != gene_col
    ]

    if (length(sample_col_candidates) > 0 && length(value_col_candidates) > 0) {
      for (sample_col in sample_col_candidates) {
        for (value_col in value_col_candidates) {
          tmp <- df %>%
            dplyr::transmute(
              gene = clean_gene_symbol(.data[[gene_col]]),
              sample_id = clean_sample_id(.data[[sample_col]]),
              value = suppressWarnings(as.numeric(.data[[value_col]]))
            ) %>%
            dplyr::filter(gene %in% required_genes, !is.na(sample_id), !is.na(value))

          if (n_distinct(tmp$gene) == length(required_genes) &&
              n_distinct(tmp$sample_id) >= expected_n_total) {
            wide <- tmp %>%
              dplyr::group_by(gene, sample_id) %>%
              dplyr::summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
              tidyr::pivot_wider(names_from = sample_id, values_from = value)

            mat <- wide %>%
              tibble::column_to_rownames("gene") %>%
              as.matrix()
            storage.mode(mat) <- "numeric"

            return(list(
              matrix = mat,
              source_file = source_file,
              format = "long",
              gene_col = gene_col,
              sample_col = sample_col,
              value_col = value_col,
              genes_present = rownames(mat)
            ))
          }
        }
      }
    }
  }

  # Wide-format support: one gene column and sample columns.
  if (length(gene_col_candidates) == 0) {
    return(NULL)
  }

  gene_counts <- vapply(gene_col_candidates, function(gc) {
    sum(clean_gene_symbol(df[[gc]]) %in% required_genes, na.rm = TRUE)
  }, numeric(1))
  gene_col <- gene_col_candidates[which.max(gene_counts)]

  non_gene_cols <- setdiff(cn, gene_col)
  numeric_counts <- vapply(df[non_gene_cols], numeric_column_score, numeric(1))
  numeric_cols <- non_gene_cols[numeric_counts >= length(required_genes)]

  if (length(numeric_cols) < expected_n_total) {
    return(NULL)
  }

  expr <- df[, c(gene_col, numeric_cols), drop = FALSE]
  names(expr)[1] <- "gene"
  expr$gene <- clean_gene_symbol(expr$gene)

  expr_required <- expr %>%
    dplyr::filter(gene %in% required_genes)

  if (n_distinct(expr_required$gene) < length(required_genes)) {
    return(NULL)
  }

  expr_required <- expr_required %>%
    dplyr::mutate(dplyr::across(all_of(numeric_cols), ~ suppressWarnings(as.numeric(.x)))) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(dplyr::across(all_of(numeric_cols), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

  mat <- expr_required %>%
    tibble::column_to_rownames("gene") %>%
    as.matrix()
  storage.mode(mat) <- "numeric"

  return(list(
    matrix = mat,
    source_file = source_file,
    format = "wide",
    gene_col = gene_col,
    sample_col = NA_character_,
    value_col = NA_character_,
    genes_present = rownames(mat)
  ))
}

safe_unpack_archives <- function(dir_path) {
  archive_files <- list.files(
    dir_path,
    pattern = "\\.(tar|tar\\.gz|tgz|zip)$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(archive_files) == 0) return(invisible(NULL))

  for (f in archive_files) {
    message("Unpacking archive: ", f)
    target_dir <- file.path(dirname(f), tools::file_path_sans_ext(basename(f)))
    dir.create(target_dir, showWarnings = FALSE, recursive = TRUE)

    tryCatch({
      if (stringr::str_detect(stringr::str_to_lower(f), "\\.zip$")) {
        utils::unzip(f, exdir = target_dir)
      } else {
        utils::untar(f, exdir = target_dir)
      }
    }, error = function(e) {
      warning("Could not unpack archive: ", f, " | ", conditionMessage(e))
    })
  }

  invisible(NULL)
}

find_expression_data <- function() {
  # 1. Local expression file has priority.
  if (file.exists(local_expression_file)) {
    df <- read_table_flexible(local_expression_file)
    obj <- detect_expression_matrix(df, local_expression_file)
    if (is.null(obj)) {
      stop("Local expression file exists but did not contain the required four genes in a usable matrix: ",
           local_expression_file, call. = FALSE)
    }
    return(obj)
  }

  # 2. Scan working directory and output directory for candidate files.
  local_candidates <- unique(c(
    list.files(".", pattern = "\\.(csv|tsv|txt)(\\.gz)?$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE),
    list.files(outdir, pattern = "\\.(csv|tsv|txt)(\\.gz)?$", full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  ))

  local_candidates <- local_candidates[
    stringr::str_detect(
      basename(local_candidates),
      stringr::regex("GSE235863|TPM|tpm|bulk|RNA|rna|expr|expression|gene", ignore_case = TRUE)
    )
  ]

  if (length(local_candidates) > 0) {
    for (f in local_candidates) {
      obj <- tryCatch({
        df <- read_table_flexible(f)
        detect_expression_matrix(df, f)
      }, error = function(e) NULL)

      if (!is.null(obj)) return(obj)
    }
  }

  # 3. Download processed GEO supplementary files if allowed.
  if (!allow_geo_download) {
    stop("No usable local expression matrix found and allow_geo_download is FALSE.", call. = FALSE)
  }

  if (!geoquery_available) {
    stop(
      "No usable local expression matrix found and GEOquery is unavailable.\n",
      "Install GEOquery or place a TPM matrix at: ", local_expression_file,
      call. = FALSE
    )
  }

  message("No local TPM matrix detected. Attempting GEO supplementary download for ", gse_id, ".")

  geo_dir <- file.path(outdir, "GEO_download")
  dir.create(geo_dir, showWarnings = FALSE, recursive = TRUE)

  # Filter aims to avoid downloading unrelated large h5ad files when possible.
  # If no file matches the filter, manually download the processed TPM matrix
  # from GEO and save it as local_expression_file.
  tryCatch({
    GEOquery::getGEOSuppFiles(
      gse_id,
      makeDirectory = TRUE,
      baseDir = geo_dir,
      fetch_files = TRUE,
      filter_regex = "(?i)(bulk|TPM|tpm|rna|RNA|gene|expr|expression|count|counts|txt|csv)"
    )
  }, error = function(e) {
    warning("GEO supplementary download produced an error: ", conditionMessage(e))
  })

  safe_unpack_archives(geo_dir)

  downloaded_candidates <- list.files(
    geo_dir,
    pattern = "\\.(csv|tsv|txt)(\\.gz)?$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  downloaded_candidates <- downloaded_candidates[
    stringr::str_detect(
      basename(downloaded_candidates),
      stringr::regex("GSE235863|TPM|tpm|bulk|RNA|rna|expr|expression|gene|count", ignore_case = TRUE)
    )
  ]

  if (length(downloaded_candidates) == 0) {
    stop(
      "GEO download did not yield a candidate text/CSV TPM matrix.\n",
      "Manually download the processed GSE235863 bulk RNA-seq TPM matrix from GEO and save it as: ",
      local_expression_file,
      call. = FALSE
    )
  }

  for (f in downloaded_candidates) {
    obj <- tryCatch({
      df <- read_table_flexible(f)
      detect_expression_matrix(df, f)
    }, error = function(e) NULL)

    if (!is.null(obj)) return(obj)
  }

  stop(
    "No downloaded candidate contained all four required genes: ",
    paste(required_genes, collapse = ", "),
    call. = FALSE
  )
}

fetch_geo_pheno <- function() {
  pheno_local <- NULL

  if (file.exists(local_metadata_file)) {
    message("Reading local metadata file: ", local_metadata_file)
    pheno_local <- read_table_flexible(local_metadata_file)

    names(pheno_local) <- stringr::str_replace_all(names(pheno_local), "\\s+", "_")

    if (!"sample_id" %in% names(pheno_local)) {
      stop("Local metadata file must contain a sample_id column.", call. = FALSE)
    }
    if (!"response" %in% names(pheno_local)) {
      stop("Local metadata file must contain a response column.", call. = FALSE)
    }

    pheno_local <- pheno_local %>%
      dplyr::mutate(
        sample_id = clean_sample_id(sample_id),
        patient = if ("patient" %in% names(.)) toupper(as.character(patient)) else extract_patient_id(sample_id),
        response_raw = as.character(response),
        response_class = standardize_response(response_raw, response_field_only = TRUE),
        metadata_text = paste(dplyr::across(dplyr::everything()), collapse = " | ")
      ) %>%
      dplyr::select(sample_id, patient, response_raw, response_class, metadata_text)

    return(pheno_local)
  }

  if (!geoquery_available || !biobase_available) {
    warning("GEOquery/Biobase unavailable and no local metadata file found.")
    return(NULL)
  }

  message("Fetching GEO sample metadata for ", gse_id, ".")

  pheno <- tryCatch({
    gset <- GEOquery::getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)
    if (is.list(gset)) {
      Biobase::pData(gset[[1]])
    } else {
      Biobase::pData(gset)
    }
  }, error = function(e) {
    warning("Could not retrieve GEO series-matrix metadata: ", conditionMessage(e))
    NULL
  })

  if (is.null(pheno) || nrow(pheno) == 0) return(NULL)

  pheno <- as.data.frame(pheno, check.names = FALSE)
  pheno$sample_id <- if ("geo_accession" %in% names(pheno)) pheno$geo_accession else rownames(pheno)

  pheno$metadata_text <- apply(pheno, 1, function(z) {
    paste(as.character(z), collapse = " | ")
  })

  response_field <- extract_response_field(pheno$metadata_text)
  response_from_field <- standardize_response(response_field, response_field_only = TRUE)
  response_from_text  <- standardize_response(pheno$metadata_text, response_field_only = FALSE)

  response_class <- ifelse(!is.na(response_from_field), response_from_field, response_from_text)

  pheno_out <- pheno %>%
    dplyr::mutate(
      sample_id = clean_sample_id(sample_id),
      patient = extract_patient_id(metadata_text),
      response_raw = ifelse(!is.na(response_field), response_field, metadata_text),
      response_class = response_class
    ) %>%
    dplyr::select(sample_id, patient, response_raw, response_class, metadata_text)

  pheno_out
}

match_metadata_to_expression <- function(sample_ids, pheno) {
  sample_ids <- clean_sample_id(sample_ids)

  sample_tbl <- tibble::tibble(
    sample_id = sample_ids,
    patient_from_sample = extract_patient_id(sample_ids),
    response_from_sample = standardize_response(sample_ids, response_field_only = FALSE)
  )

  if (is.null(pheno) || nrow(pheno) == 0) {
    return(sample_tbl %>%
             dplyr::transmute(
               sample_id,
               patient = patient_from_sample,
               response_raw = sample_id,
               response_class = response_from_sample,
               metadata_text = sample_id
             ))
  }

  pheno <- pheno %>%
    dplyr::mutate(
      sample_id = clean_sample_id(sample_id),
      patient = toupper(as.character(patient))
    )

  exact <- sample_tbl %>%
    dplyr::left_join(pheno, by = "sample_id")

  # Fill missing rows by patient ID if possible.
  need_patient <- is.na(exact$response_class) & !is.na(exact$patient_from_sample)

  if (any(need_patient)) {
    pheno_patient <- pheno %>%
      dplyr::filter(!is.na(patient), !is.na(response_class)) %>%
      dplyr::group_by(patient) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::rename(patient_from_sample = patient)

    by_patient <- sample_tbl %>%
      dplyr::left_join(pheno_patient, by = "patient_from_sample")

    fill_cols <- c("patient", "response_raw", "response_class", "metadata_text")
    for (cc in fill_cols) {
      exact[[cc]][need_patient] <- by_patient[[cc]][need_patient]
    }
  }

  # Last fallback: response embedded in expression sample ID.
  exact <- exact %>%
    dplyr::mutate(
      patient = ifelse(is.na(patient), patient_from_sample, patient),
      response_raw = ifelse(is.na(response_raw), sample_id, response_raw),
      response_class = ifelse(is.na(response_class), response_from_sample, response_class),
      metadata_text = ifelse(is.na(metadata_text), sample_id, metadata_text)
    ) %>%
    dplyr::select(sample_id, patient, response_raw, response_class, metadata_text)

  exact
}

make_qc_row <- function(check, observed, expected, pass, tolerance = NA_character_) {
  tibble::tibble(
    check = check,
    observed = as.character(observed),
    expected = as.character(expected),
    tolerance = as.character(tolerance),
    pass = as.logical(pass)
  )
}

# ----------------------------
# Load expression and metadata
# ----------------------------
expr_obj <- find_expression_data()
expr_mat <- expr_obj$matrix

# Keep only required genes in coefficient order.
expr_mat <- expr_mat[required_genes, , drop = FALSE]
colnames(expr_mat) <- clean_sample_id(colnames(expr_mat))

pheno <- fetch_geo_pheno()
sample_meta <- match_metadata_to_expression(colnames(expr_mat), pheno)

# Build sample-level source data.
expr_long <- as.data.frame(t(expr_mat), check.names = FALSE) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::mutate(sample_id = clean_sample_id(sample_id))

source_df <- expr_long %>%
  dplyr::left_join(sample_meta, by = "sample_id") %>%
  dplyr::filter(!is.na(response_class)) %>%
  dplyr::mutate(
    response_class = factor(response_class, levels = c("Responder", "NonResponder")),
    non_responder = ifelse(response_class == "NonResponder", 1L, 0L)
  )

# Retain samples with complete four-gene expression.
for (g in required_genes) {
  source_df[[g]] <- suppressWarnings(as.numeric(source_df[[g]]))
}

source_df <- source_df %>%
  dplyr::filter(stats::complete.cases(dplyr::across(all_of(required_genes))))

if (nrow(source_df) == 0) {
  stop("No samples remained after response annotation and complete four-gene filtering.", call. = FALSE)
}

# Apply locked expression transform.
if (input_expression_already_log2_tpm) {
  for (g in required_genes) {
    source_df[[paste0(g, "_log2TPM1")]] <- source_df[[g]]
  }
  expression_transform_note <- "Input treated as already log2(TPM + 1)."
} else {
  if (any(as.matrix(source_df[, required_genes, drop = FALSE]) < 0, na.rm = TRUE)) {
    stop("Negative expression values detected. Locked analysis expects TPM >= 0 before log2(TPM + 1).",
         call. = FALSE)
  }
  for (g in required_genes) {
    source_df[[paste0(g, "_log2TPM1")]] <- log2(source_df[[g]] + 1)
  }
  expression_transform_note <- "log2(TPM + 1) was applied to raw TPM values."
}

score_terms <- paste0(required_genes, "_log2TPM1")

source_df <- source_df %>%
  dplyr::mutate(
    four_gene_score =
      score_coefficients["TPI1"]   * .data[["TPI1_log2TPM1"]] +
      score_coefficients["ENO1"]   * .data[["ENO1_log2TPM1"]] +
      score_coefficients["LDHA"]   * .data[["LDHA_log2TPM1"]] +
      score_coefficients["SLC2A1"] * .data[["SLC2A1_log2TPM1"]]
  ) %>%
  dplyr::mutate(
    four_gene_score_z = as.numeric(scale(four_gene_score)),
    median_score_cut = median(four_gene_score, na.rm = TRUE),
    # Locked Fisher result uses median sample included in High group.
    score_group_median_in_high = ifelse(four_gene_score >= median_score_cut, "High", "Low"),
    score_group_median_in_high = factor(score_group_median_in_high, levels = c("Low", "High"))
  ) %>%
  dplyr::arrange(response_class, dplyr::desc(four_gene_score))

# ----------------------------
# Statistics
# ----------------------------
roc_obj <- pROC::roc(
  response = source_df$non_responder,
  predictor = source_df$four_gene_score,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

auc_value <- as.numeric(pROC::auc(roc_obj))
auc_ci <- as.numeric(pROC::ci.auc(roc_obj, method = "delong"))

roc_coords <- pROC::coords(
  roc_obj,
  x = "all",
  ret = c("threshold", "specificity", "sensitivity"),
  transpose = FALSE
) %>%
  as.data.frame() %>%
  dplyr::mutate(
    fpr = 1 - specificity,
    tpr = sensitivity,
    auc = auc_value,
    ci_low = auc_ci[1],
    ci_high = auc_ci[3],
    positive_class = "NonResponder"
  ) %>%
  dplyr::arrange(fpr, tpr)

strat_table <- source_df %>%
  dplyr::count(score_group_median_in_high, response_class, name = "n") %>%
  tidyr::complete(
    score_group_median_in_high = factor(c("Low", "High"), levels = c("Low", "High")),
    response_class = factor(c("Responder", "NonResponder"), levels = c("Responder", "NonResponder")),
    fill = list(n = 0L)
  )

nr_high <- strat_table %>%
  dplyr::filter(score_group_median_in_high == "High", response_class == "NonResponder") %>%
  dplyr::pull(n)
r_high <- strat_table %>%
  dplyr::filter(score_group_median_in_high == "High", response_class == "Responder") %>%
  dplyr::pull(n)
nr_low <- strat_table %>%
  dplyr::filter(score_group_median_in_high == "Low", response_class == "NonResponder") %>%
  dplyr::pull(n)
r_low <- strat_table %>%
  dplyr::filter(score_group_median_in_high == "Low", response_class == "Responder") %>%
  dplyr::pull(n)

fisher_matrix <- matrix(
  c(nr_high, r_high, nr_low, r_low),
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    score_group = c("High", "Low"),
    response = c("NonResponder", "Responder")
  )
)

fisher_p <- stats::fisher.test(fisher_matrix, alternative = "greater")$p.value

strat_table_out <- as.data.frame(as.table(fisher_matrix)) %>%
  dplyr::rename(n = Freq) %>%
  dplyr::mutate(
    median_rule = "High if score >= median; median sample included in High group",
    fisher_alternative = "greater enrichment of non-responders in High group",
    fisher_p = fisher_p
  )

# ----------------------------
# Write source CSVs before QC/figures
# ----------------------------
manifest <- tibble::tibble(
  item = c(
    "analysis_script",
    "gse_id",
    "expression_source_file",
    "expression_source_format",
    "expression_transform",
    "score_formula",
    "positive_class",
    "interpretation"
  ),
  value = c(
    "12_drug_repurposing.R",
    gse_id,
    expr_obj$source_file,
    expr_obj$format,
    expression_transform_note,
    "0.3041908*TPI1 + 0.9639654*ENO1 + 1.3404374*LDHA + 0.2424239*SLC2A1",
    "NonResponder",
    "Hypothesis-generating exploratory association only; not a treatment-selection assay."
  )
)

source_csv <- file.path(outdir, "GSE235863_four_gene_score_source_data.csv")
roc_source_csv <- file.path(outdir, "GSE235863_ROC_curve_source_data.csv")
strat_source_csv <- file.path(outdir, "GSE235863_median_stratification_source_data.csv")
manifest_csv <- file.path(outdir, "GSE235863_input_manifest.csv")

readr::write_csv(source_df, source_csv)
readr::write_csv(roc_coords, roc_source_csv)
readr::write_csv(strat_table_out, strat_source_csv)
readr::write_csv(manifest, manifest_csv)

message("Wrote source CSV: ", source_csv)
message("Wrote ROC source CSV: ", roc_source_csv)
message("Wrote stratification source CSV: ", strat_source_csv)
message("Wrote manifest CSV: ", manifest_csv)

# ----------------------------
# Mandatory QC
# ----------------------------
n_total <- nrow(source_df)
n_resp <- sum(source_df$response_class == "Responder")
n_nr <- sum(source_df$response_class == "NonResponder")
genes_present <- all(required_genes %in% rownames(expr_mat))
all_nr_high <- all(
  source_df$score_group_median_in_high[source_df$response_class == "NonResponder"] == "High"
)

qc <- dplyr::bind_rows(
  make_qc_row("GSE accession", gse_id, "GSE235863", identical(gse_id, "GSE235863")),
  make_qc_row("Required genes present", paste(required_genes, collapse = ";"),
              "TPI1;ENO1;LDHA;SLC2A1", genes_present),
  make_qc_row("Total evaluable samples", n_total, expected_n_total, n_total == expected_n_total),
  make_qc_row("Responders", n_resp, expected_n_responder, n_resp == expected_n_responder),
  make_qc_row("Non-responders", n_nr, expected_n_nonresponder, n_nr == expected_n_nonresponder),
  make_qc_row("Positive ROC class", "NonResponder", "NonResponder", TRUE),
  make_qc_row("Expression transform", expression_transform_note,
              "log2(TPM + 1) used for score", TRUE),
  make_qc_row("ROC AUC", sprintf("%.6f", auc_value), expected_auc,
              abs(auc_value - expected_auc) <= auc_tolerance, auc_tolerance),
  make_qc_row("ROC 95% CI lower", sprintf("%.6f", auc_ci[1]), expected_ci_low,
              abs(auc_ci[1] - expected_ci_low) <= ci_tolerance, ci_tolerance),
  make_qc_row("ROC 95% CI upper", sprintf("%.6f", auc_ci[3]), expected_ci_high,
              abs(auc_ci[3] - expected_ci_high) <= ci_tolerance, ci_tolerance),
  make_qc_row("All non-responders in median-high group", all_nr_high, TRUE, all_nr_high),
  make_qc_row("Fisher p, median sample included in High", sprintf("%.6f", fisher_p),
              expected_fisher_p, abs(fisher_p - expected_fisher_p) <= fisher_tolerance,
              fisher_tolerance),
  make_qc_row("Interpretation label",
              "hypothesis-generating exploratory only",
              "hypothesis-generating exploratory only",
              TRUE)
)

qc_status <- ifelse(all(qc$pass), "PASS", "FAIL")
qc <- qc %>%
  dplyr::mutate(overall_qc = qc_status)

qc_csv <- file.path(outdir, "GSE235863_nonresponse_QC_check.csv")
readr::write_csv(qc, qc_csv)
message("Wrote QC CSV: ", qc_csv)

print(qc)

if (!all(qc$pass)) {
  stop(
    "Mandatory QC failed. Source CSV and QC CSV have been written, but PDF/PNG figures were not generated.\n",
    "Inspect: ", qc_csv,
    call. = FALSE
  )
}

# ----------------------------
# Plot only after QC PASS
# ----------------------------
auc_label <- paste0(
  "AUC = ", sprintf("%.3f", auc_value), "\n",
  "95% CI = ", sprintf("%.3f", auc_ci[1]), "-", sprintf("%.3f", auc_ci[3]), "\n",
  "Positive class: non-response\n",
  "n = ", n_total, " (R = ", n_resp, ", NR = ", n_nr, ")"
)

roc_plot <- ggplot2::ggplot(roc_coords, ggplot2::aes(x = fpr, y = tpr)) +
  ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.4) +
  ggplot2::geom_step(linewidth = 0.9, direction = "hv") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  ggplot2::annotate(
    "text",
    x = 0.62,
    y = 0.18,
    label = auc_label,
    hjust = 0,
    vjust = 0,
    size = 3.4
  ) +
  ggplot2::labs(
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(color = "black"),
    plot.margin = ggplot2::margin(8, 10, 8, 8)
  )

fig_pdf <- file.path(outdir, "Fig7B_GSE235863_four_gene_nonresponse_ROC.pdf")
fig_png <- file.path(outdir, "Fig7B_GSE235863_four_gene_nonresponse_ROC.png")

ggplot2::ggsave(fig_pdf, roc_plot, width = 4.2, height = 4.0, useDingbats = FALSE)
ggplot2::ggsave(fig_png, roc_plot, width = 4.2, height = 4.0, dpi = 600)

# A compact score-distribution panel is saved as a supplementary QC visualization,
# not as an additional manuscript claim.
score_plot <- ggplot2::ggplot(
  source_df,
  ggplot2::aes(x = response_class, y = four_gene_score)
) +
  ggplot2::geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.4) +
  ggplot2::geom_point(
    ggplot2::aes(shape = score_group_median_in_high),
    size = 2.4,
    position = ggplot2::position_jitter(width = 0.08, height = 0)
  ) +
  ggplot2::geom_hline(yintercept = median(source_df$four_gene_score), linetype = "dashed", linewidth = 0.4) +
  ggplot2::labs(
    x = NULL,
    y = "Four-gene glycolysis score",
    shape = "Median group"
  ) +
  ggplot2::theme_classic(base_size = 12) +
  ggplot2::theme(
    axis.text = ggplot2::element_text(color = "black"),
    axis.title = ggplot2::element_text(color = "black"),
    legend.position = "right"
  )

score_pdf <- file.path(outdir, "GSE235863_four_gene_score_by_response_QC_plot.pdf")
score_png <- file.path(outdir, "GSE235863_four_gene_score_by_response_QC_plot.png")

ggplot2::ggsave(score_pdf, score_plot, width = 4.4, height = 3.6, useDingbats = FALSE)
ggplot2::ggsave(score_png, score_plot, width = 4.4, height = 3.6, dpi = 600)

message("QC PASS. Figure files written:")
message("  ", fig_pdf)
message("  ", fig_png)
message("  ", score_pdf)
message("  ", score_png)

message("Completed 12_drug_repurposing.R replacement analysis.")
