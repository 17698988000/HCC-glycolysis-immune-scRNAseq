# Helper configuration for revision figure restoration scripts
# Usage: source("00_config_paths.R")

project_dir <- Sys.getenv("PROJECT_DIR", unset = "D:/scRNA_project")
if (!dir.exists(project_dir)) {
  stop("PROJECT_DIR does not exist: ", project_dir,
       "\nSet Sys.setenv(PROJECT_DIR='path/to/scRNA_project') before running.")
}
setwd(project_dir)

ensure_packages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("Missing required R packages: ", paste(missing, collapse = ", "),
         "\nInstall them before running this script.")
  }
  invisible(TRUE)
}

fmt_p <- function(p) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "p < 0.001", paste0("p = ", sprintf("%.3f", p))))
}
