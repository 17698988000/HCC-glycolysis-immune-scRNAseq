# Run selected locked figure restoration scripts.
#
# Usage:
#   From repository root:
#     source("revision_figure_restore/00_run_selected_locked_figures.R")
#
#   From inside revision_figure_restore/:
#     source("00_run_selected_locked_figures.R")
#
#   Or explicitly:
#     Sys.setenv(PROJECT_DIR = "/path/to/HCC-glycolysis-immune-scRNAseq")
#     source("revision_figure_restore/00_run_selected_locked_figures.R")
#
# This runner intentionally uses full script paths and chdir = TRUE so that each
# figure script can still source("00_config_paths.R") from its own directory.

# -----------------------------------------------------------------------------
# Locate and source shared config
# -----------------------------------------------------------------------------

config_candidates <- c(
  file.path("revision_figure_restore", "00_config_paths.R"),
  "00_config_paths.R"
)

config_file <- config_candidates[file.exists(config_candidates)][1]

if (is.na(config_file)) {
  stop(
    "Could not find 00_config_paths.R. Run from repository root, ",
    "run from revision_figure_restore/, or set PROJECT_DIR before running."
  )
}

source(config_file, local = FALSE, chdir = TRUE)

# 00_config_paths.R should define revision_figure_restore_dir.
if (!exists("revision_figure_restore_dir")) {
  stop("00_config_paths.R did not define revision_figure_restore_dir.")
}

# -----------------------------------------------------------------------------
# Selected locked figure restoration scripts
# -----------------------------------------------------------------------------

scripts <- c(
  "plot_Fig4E_SC_ENO1_SPP1_correlation.R",
  "plot_Fig5B_ENO1_KM_optimal.R",
  "plot_Fig5C_ENO1_AJCC_stage.R",
  "plot_Fig5D_ENO1_Cox_forest.R",
  "plot_Fig6A_four_gene_score_KM.R",
  "plot_Fig6B_four_gene_score_Cox_forest.R",
  "plot_Fig9B_spatial_direction_consistency.R",
  "plot_FigS20_RCTD_celltype_composition.R",
  "plot_FigS24_calibration_TCGA_GSE14520.R",
  "plot_FigS25_GCK_sensitivity.R"
)

script_paths <- file.path(revision_figure_restore_dir, scripts)
missing_scripts <- scripts[!file.exists(script_paths)]

if (length(missing_scripts) > 0) {
  stop(
    "Missing figure restoration script(s): ",
    paste(missing_scripts, collapse = ", "),
    "\nExpected under: ",
    revision_figure_restore_dir
  )
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------

run_log <- data.frame(
  script = scripts,
  status = NA_character_,
  start_time = as.character(NA),
  end_time = as.character(NA),
  error = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(script_paths)) {
  s <- scripts[i]
  s_path <- script_paths[i]

  message("\n===== Running ", s, " =====")
  run_log$start_time[i] <- as.character(Sys.time())

  result <- tryCatch(
    {
      source(s_path, local = new.env(parent = globalenv()), chdir = TRUE)
      list(status = "COMPLETED", error = NA_character_)
    },
    error = function(e) {
      list(status = "FAILED", error = conditionMessage(e))
    }
  )

  run_log$status[i] <- result$status
  run_log$error[i] <- result$error
  run_log$end_time[i] <- as.character(Sys.time())

  if (identical(result$status, "FAILED")) {
    message("FAILED: ", s)
    message("Error: ", result$error)
    break
  }
}

# Save runner log if config created revision_results_dir; otherwise save locally.
log_dir <- if (exists("revision_results_dir") && dir.exists(revision_results_dir)) {
  revision_results_dir
} else {
  getwd()
}

log_file <- file.path(log_dir, "00_run_selected_locked_figures_log.csv")
utils::write.csv(run_log, log_file, row.names = FALSE, quote = TRUE)
message("\nRun log written to: ", normalizePath(log_file, winslash = "/", mustWork = FALSE))

if (any(run_log$status == "FAILED", na.rm = TRUE)) {
  stop("At least one figure restoration script failed. See run log: ", log_file)
}

message("\nAll selected locked figure restoration scripts completed.")
