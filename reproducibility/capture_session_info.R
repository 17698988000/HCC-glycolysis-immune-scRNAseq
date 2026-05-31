#!/usr/bin/env Rscript

# Capture the current R session and selected package versions.
# Run from the repository root:
#   Rscript reproducibility/capture_session_info.R

outdir <- "reproducibility"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "current_author_sessionInfo.txt"))

packages <- c(
  "R", "Seurat", "SeuratObject", "harmony", "AUCell", "CellChat", "infercnv",
  "glmnet", "survival", "survminer", "timeROC", "rms", "estimate", "spacexr",
  "nichenetr", "liana", "SingleCellExperiment", "SummarizedExperiment",
  "dplyr", "tidyr", "readr", "tibble", "ggplot2", "Matrix"
)
installed <- installed.packages()
versions <- data.frame(
  package = packages,
  version = vapply(packages, function(pkg) {
    if (pkg == "R") return(as.character(getRversion()))
    if (pkg %in% rownames(installed)) as.character(installed[pkg, "Version"]) else "NOT_INSTALLED_IN_CURRENT_SESSION"
  }, character(1)),
  stringsAsFactors = FALSE
)
utils::write.csv(versions, file.path(outdir, "current_author_package_versions.csv"), row.names = FALSE)
message("Environment snapshot written to reproducibility/.")
