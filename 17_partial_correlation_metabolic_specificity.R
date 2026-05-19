############################################################
## Supplementary Figure S21
## Metabolic specificity analysis:
## Partial Spearman correlations of immunosuppressive ligands
## with glycolysis versus OXPHOS activity
##
## Required input:
## - seurat_final.rds
## - metadata column: Glycolysis_AUC
## - metadata column: OXPHOS_AUC or OxPhos_AUC or OXPHOS_score
##
## Output:
## - FigS21_partial_correlation_metabolic_specificity.pdf
## - FigS21_partial_correlation_metabolic_specificity.png
## - FigS21_partial_correlation_metabolic_specificity_source.csv
############################################################

rm(list = ls())

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(tibble)
})

setwd("D:/scRNA_project")

outdir <- "figures_final"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Load object
############################################################

obj_file <- "seurat_final.rds"

if (!file.exists(obj_file)) {
  stop("Cannot find seurat_final.rds in D:/scRNA_project")
}

obj <- readRDS(obj_file)
meta <- obj@meta.data

message("Loaded: ", obj_file)
message("Total cells: ", ncol(obj))

############################################################
## 2. Required metadata columns
############################################################

patient_col <- "patient"
celltype_col <- "cell_type"
site_col <- "site"
gly_col <- "Glycolysis_AUC"

oxphos_candidates <- c(
  "OXPHOS_AUC",
  "OxPhos_AUC",
  "OXPHOS_score",
  "OxPhos_score",
  "oxidative_phosphorylation_AUC",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION_AUC"
)

ox_col <- oxphos_candidates[oxphos_candidates %in% colnames(meta)][1]

if (!patient_col %in% colnames(meta)) stop("Missing metadata column: patient")
if (!celltype_col %in% colnames(meta)) stop("Missing metadata column: cell_type")
if (!site_col %in% colnames(meta)) stop("Missing metadata column: site")
if (!gly_col %in% colnames(meta)) stop("Missing metadata column: Glycolysis_AUC")

if (is.na(ox_col)) {
  stop(
    paste0(
      "No OXPHOS score column found. Expected one of: ",
      paste(oxphos_candidates, collapse = ", "),
      "\nPlease run 14_OXPHOS_metabolic_specificity.R first, or add the OXPHOS AUCell score to obj@meta.data."
    )
  )
}

message("Using glycolysis column: ", gly_col)
message("Using OXPHOS column: ", ox_col)

############################################################
## 3. Select tumor-derived hepatocytes
############################################################

celltype_vec <- as.character(meta[[celltype_col]])
site_vec <- as.character(meta[[site_col]])

is_tumor_site <- site_vec == "Tumor"
is_tumor_hep <- grepl("hep|tumor", celltype_vec, ignore.case = TRUE)

tumor_hep_cells <- rownames(meta)[is_tumor_site & is_tumor_hep]

message("Selected tumor-derived hepatocytes: ", length(tumor_hep_cells))
print(table(meta[tumor_hep_cells, patient_col]))

if (length(tumor_hep_cells) != 15391) {
  warning(
    paste0(
      "Selected cells = ", length(tumor_hep_cells),
      ", not 15,391. Please verify cell_type/site filtering."
    )
  )
}

############################################################
## 4. Extract ligand expression
############################################################

ligands <- c("SPP1", "MIF")

assay_use <- "RNA"
if (!assay_use %in% Assays(obj)) assay_use <- DefaultAssay(obj)

expr_mat <- tryCatch(
  GetAssayData(obj, assay = assay_use, layer = "data"),
  error = function(e) GetAssayData(obj, assay = assay_use, slot = "data")
)

missing_ligands <- setdiff(ligands, rownames(expr_mat))
if (length(missing_ligands) > 0) {
  stop("Missing ligand genes in expression matrix: ", paste(missing_ligands, collapse = ", "))
}

expr_sub <- as.matrix(expr_mat[ligands, tumor_hep_cells, drop = FALSE])

dat <- tibble(
  cell = tumor_hep_cells,
  patient = as.character(meta[tumor_hep_cells, patient_col]),
  glycolysis = as.numeric(meta[tumor_hep_cells, gly_col]),
  oxphos = as.numeric(meta[tumor_hep_cells, ox_col]),
  SPP1 = as.numeric(expr_sub["SPP1", ]),
  MIF = as.numeric(expr_sub["MIF", ])
)

write_csv(
  dat,
  file.path(outdir, "FigS21_partial_correlation_input_cells.csv")
)

############################################################
## 5. Partial Spearman correlation helper
##
## Partial Spearman is computed by:
## 1. ranking x, y, z;
## 2. residualizing ranked x and ranked y against ranked z;
## 3. correlating residuals.
############################################################

partial_spearman <- function(x, y, z) {
  ok <- complete.cases(x, y, z)
  x <- x[ok]
  y <- y[ok]
  z <- z[ok]
  
  n <- length(x)
  
  if (n < 10) {
    return(tibble(
      n = n,
      rho = NA_real_,
      p_value = NA_real_
    ))
  }
  
  xr <- rank(x, ties.method = "average")
  yr <- rank(y, ties.method = "average")
  zr <- rank(z, ties.method = "average")
  
  rx <- residuals(lm(xr ~ zr))
  ry <- residuals(lm(yr ~ zr))
  
  rho <- cor(rx, ry, method = "pearson")
  
  ## Approximate test for partial correlation with one covariate
  df <- n - 3
  tval <- rho * sqrt(df / (1 - rho^2))
  pval <- 2 * pt(abs(tval), df = df, lower.tail = FALSE)
  
  tibble(
    n = n,
    rho = rho,
    p_value = pval
  )
}

############################################################
## 6. Run S21 partial correlation analysis
############################################################

res <- bind_rows(
  partial_spearman(dat$SPP1, dat$glycolysis, dat$oxphos) %>%
    mutate(
      ligand = "SPP1",
      comparison = "Glycolysis adjusted for OXPHOS"
    ),
  partial_spearman(dat$SPP1, dat$oxphos, dat$glycolysis) %>%
    mutate(
      ligand = "SPP1",
      comparison = "OXPHOS adjusted for glycolysis"
    ),
  partial_spearman(dat$MIF, dat$glycolysis, dat$oxphos) %>%
    mutate(
      ligand = "MIF",
      comparison = "Glycolysis adjusted for OXPHOS"
    ),
  partial_spearman(dat$MIF, dat$oxphos, dat$glycolysis) %>%
    mutate(
      ligand = "MIF",
      comparison = "OXPHOS adjusted for glycolysis"
    )
) %>%
  mutate(
    ligand = factor(ligand, levels = c("SPP1", "MIF")),
    comparison = factor(
      comparison,
      levels = c(
        "Glycolysis adjusted for OXPHOS",
        "OXPHOS adjusted for glycolysis"
      )
    ),
    FDR = p.adjust(p_value, method = "BH"),
    label = sprintf("%.2f", rho)
  )

print(res)

write_csv(
  res,
  file.path(outdir, "FigS21_partial_correlation_metabolic_specificity_source.csv")
)

############################################################
## 7. Plot S21
############################################################

p <- ggplot(
  res,
  aes(x = ligand, y = rho, fill = comparison)
) +
  geom_hline(yintercept = 0, linewidth = 0.4, color = "grey50") +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    color = "black",
    linewidth = 0.25
  ) +
  geom_text(
    aes(label = label),
    position = position_dodge(width = 0.75),
    vjust = ifelse(res$rho >= 0, -0.45, 1.25),
    size = 4
  ) +
  coord_cartesian(ylim = c(
    min(-0.15, min(res$rho, na.rm = TRUE) - 0.08),
    max(0.75, max(res$rho, na.rm = TRUE) + 0.08)
  )) +
  labs(
    title = "Metabolic specificity of ligand associations",
    subtitle = "Partial Spearman correlations in tumor-derived hepatocytes",
    x = NULL,
    y = "Partial Spearman rho",
    fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "bottom"
  )

ggsave(
  filename = file.path(outdir, "FigS21_partial_correlation_metabolic_specificity.pdf"),
  plot = p,
  width = 6.8,
  height = 4.6,
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "FigS21_partial_correlation_metabolic_specificity.png"),
  plot = p,
  width = 6.8,
  height = 4.6,
  dpi = 600
)

print(p)

message("Done. Outputs written to: ", normalizePath(outdir))
