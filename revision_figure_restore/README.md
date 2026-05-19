
# Final figure restoration scripts

This folder contains final submission-ready clean vector figure scripts generated after figure QC.

The original analysis scripts in the repository remain the primary workflow. The scripts in this folder regenerate selected final manuscript figures after figure-level quality control, renumbering, and vector-output restoration.

## Included scripts

- `plot_Fig4E_SC_ENO1_SPP1_correlation.R`: single-cell ENO1-SPP1 correlation.
- `plot_Fig5B_ENO1_KM_optimal.R`: ENO1 Kaplan-Meier survival curve.
- `plot_Fig5C_ENO1_AJCC_stage.R`: ENO1 expression across AJCC pathologic stage.
- `plot_Fig5D_ENO1_Cox_forest.R`: ENO1 multivariate Cox forest plot.
- `plot_Fig6A_four_gene_score_KM.R`: four-gene score Kaplan-Meier survival curve.
- `plot_Fig6B_four_gene_score_Cox_forest.R`: four-gene score multivariate Cox forest plot.
- `plot_Fig9B_spatial_direction_consistency.R`: spatial GlycoHigh versus GlycoLow direction consistency. This is Figure 9B, not Figure 8.
- `plot_FigS20_RCTD_celltype_composition.R`: RCTD-estimated spatial cell-type composition.
- `plot_FigS24_calibration_TCGA_GSE14520.R`: TCGA-LIHC and GSE14520 3-year OS calibration.
- `plot_FigS25_GCK_sensitivity.R`: GCK sensitivity analysis.

## Notes

- Figure 8 is reserved for benchmarking / model comparison and is not regenerated in this folder yet.
- Figure 9 is the spatial transcriptomics figure.
- Supplementary Figure S20 is RCTD cell-type composition.
- Supplementary Figure S22 is NicheNet ligand activity.
- Supplementary Figure S23 is patient-level ligand mean-expression robustness.
- Supplementary Figure S24 is calibration.
- Supplementary Figure S25 is GCK sensitivity.
- Tentative figures, including the current FigS15 glycolysis gradient candidate, are not included here until their data scope is finalized.
