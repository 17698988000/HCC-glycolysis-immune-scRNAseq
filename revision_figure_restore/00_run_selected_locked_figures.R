# Run selected locked figure restoration scripts.
# Run from repository root or set PROJECT_DIR first.

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

for (s in scripts) {
  message("\n===== Running ", s, " =====")
  source(s, local = new.env(parent = globalenv()))
}
