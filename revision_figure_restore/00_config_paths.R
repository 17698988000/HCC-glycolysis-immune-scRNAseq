# Helper configuration for revision figure restoration scripts
#
# Usage options:
#   1. From project root:
#        source("revision_figure_restore/00_config_paths.R")
#
#   2. From inside revision_figure_restore/:
#        source("00_config_paths.R")
#
#   3. Explicitly set project directory before sourcing:
#        Sys.setenv(PROJECT_DIR = "/path/to/HCC-glycolysis-immune-scRNAseq")
#        source("revision_figure_restore/00_config_paths.R")
#
# This file intentionally has no machine-specific default such as D:/scRNA_project.

# -----------------------------------------------------------------------------
# Project directory detection
# -----------------------------------------------------------------------------

detect_project_dir <- function() {
  env_dir <- Sys.getenv("PROJECT_DIR", unset = NA_character_)

  if (!is.na(env_dir) && nzchar(env_dir)) {
    return(normalizePath(env_dir, winslash = "/", mustWork = FALSE))
  }

  cwd <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

  # Case A: user is already in revision_figure_restore/
  if (basename(cwd) == "revision_figure_restore") {
    return(normalizePath(dirname(cwd), winslash = "/", mustWork = TRUE))
  }

  # Case B: user is in project root
  if (dir.exists(file.path(cwd, "revision_figure_restore"))) {
    return(cwd)
  }

  # Case C: user is one level below project root for some reason
  parent <- normalizePath(file.path(cwd, ".."), winslash = "/", mustWork = TRUE)
  if (dir.exists(file.path(parent, "revision_figure_restore"))) {
    return(parent)
  }

  stop(
    "Could not infer project directory from current working directory:\n  ",
    cwd,
    "\nRun from the repository root, run from revision_figure_restore/, ",
    "or set Sys.setenv(PROJECT_DIR='/path/to/HCC-glycolysis-immune-scRNAseq') before sourcing."
  )
}

project_dir <- detect_project_dir()

if (!dir.exists(project_dir)) {
  stop("PROJECT_DIR does not exist: ", project_dir)
}

revision_figure_restore_dir <- file.path(project_dir, "revision_figure_restore")

if (!dir.exists(revision_figure_restore_dir)) {
  stop(
    "revision_figure_restore directory not found under PROJECT_DIR:\n  ",
    revision_figure_restore_dir
  )
}

# Work from project root so all scripts use consistent relative paths.
setwd(project_dir)

message("PROJECT_DIR: ", normalizePath(project_dir, winslash = "/", mustWork = TRUE))
message("revision_figure_restore_dir: ",
        normalizePath(revision_figure_restore_dir, winslash = "/", mustWork = TRUE))

# -----------------------------------------------------------------------------
# Common output directories
# -----------------------------------------------------------------------------

results_dir <- file.path(project_dir, "results")
figure_dir <- file.path(project_dir, "figures")
revision_results_dir <- file.path(results_dir, "revision_figure_restore")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(revision_results_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing required R packages: ",
      paste(missing, collapse = ", "),
      "\nInstall them before running this script."
    )
  }
  invisible(TRUE)
}

fmt_p <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p)))
  )
}

fmt_p_sci <- function(p) {
  ifelse(
    is.na(p),
    NA_character_,
    ifelse(p < 0.001, format(p, scientific = TRUE, digits = 3),
           paste0("p = ", sprintf("%.3f", p)))
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, sprintf(paste0("%.", digits, "f"), x))
}

safe_read_csv <- function(path, ...) {
  if (!file.exists(path)) {
    stop("Required file not found: ", path)
  }
  utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

safe_save_pdf_png <- function(plot, basename_no_ext, width, height, dpi = 300) {
  ensure_packages(c("ggplot2"))

  pdf_path <- file.path(figure_dir, paste0(basename_no_ext, ".pdf"))
  png_path <- file.path(figure_dir, paste0(basename_no_ext, ".png"))

  ggplot2::ggsave(pdf_path, plot = plot, width = width, height = height)
  ggplot2::ggsave(png_path, plot = plot, width = width, height = height, dpi = dpi)

  message("Wrote: ", normalizePath(pdf_path, winslash = "/", mustWork = FALSE))
  message("Wrote: ", normalizePath(png_path, winslash = "/", mustWork = FALSE))

  invisible(c(pdf = pdf_path, png = png_path))
}
