
# HCC-glycolysis-immune-scRNAseq

Analysis code for ENO1/glycolysis-immune evasion study in HCC

## Scripts

| Script | Description |
|---|---|
| `01_QC_clustering.R` | Quality control, Harmony batch correction, and cell type annotation (Sections 2.1–2.3) |
| `02_glycolysis_scoring.R` | AUCell glycolysis scoring and GlycoHigh/GlycoLow stratification (Section 2.4) |
| `03_cellchat_analysis.R` | CellChat cell-cell communication analysis (Section 2.6) |
| `04_TCGA_validation.R` | TCGA-LIHC ENO1 survival analysis, immune infiltration, and TIDE scores (Sections 2.7–2.8) |
| `05_GSE125449_validation.R` | Cross-dataset validation in GSE125449 (Section 3.8) |
| `06_inferCNV_malignant.R` | inferCNV malignant hepatocyte identification (Section 2.5) |
| `07_LASSO_risk_score.R` | LASSO-penalized Cox regression glycolysis risk score (Section 2.10) |
| `08_glycolysis_gradient.R` | Supplementary Figure S15 glycolysis-gradient analysis using 10 equal-width AUCell-score bins and mean log-normalized expression (Section 2.11; Results Section 3.12) |
| `09_TF_activity.R` | DoRothEA/VIPER transcription factor activity inference (Section 2.12) |
| `10_GSE14520_validation.R` | Independent validation of risk score in GSE14520 (Section 2.13) |
| `11_spatial_transcriptomics.R` | Visium spatial transcriptomics and RCTD deconvolution (Section 2.14) |
| `12_drug_repurposing.R` | Drug repurposing enrichment analysis (Section 2.17) |
| `13_revision_analyses.R` | Revision-stage supplementary analyses (Figures S9–S13, S20) |
| `14_OXPHOS_metabolic_specificity.R` | OXPHOS AUCell scoring and metabolic specificity analysis (Section 2.15) |
| `15_NicheNet_analysis.R` | NicheNet ligand activity analysis (Section 2.16, Supplementary Figure S22) |
| `16_FigS23_patient_level_ligand_effects.R` | Patient-level ligand mean-expression robustness analysis for SPP1, MIF, and PTGES (Supplementary Figure S23) |
| `17_FigS21_partial_correlation_metabolic_specificity.R` | Supplementary Figure S21 partial Spearman correlation plot for glycolysis versus OXPHOS specificity |

## Final figure restoration scripts

Final submission-ready clean vector figures generated after figure QC are provided in:

`revision_figure_restore/`

These scripts reproduce selected final manuscript figures after figure-level quality control, renumbering, and vector-output restoration. Original analysis scripts (`01_*.R` to `17_*.R`) remain the primary analysis workflow.

Corrected final figure mapping:

- Figure 8: benchmarking / model comparison.
- Supplementary Figure S12: per-patient ENO1-glycolysis correlation in tumor-derived hepatocytes.
- Supplementary Figure S15: glycolysis-gradient expression of ENO1, LDHA, SPP1, MIF, and PTGES in tumor-derived hepatocytes, reproduced by `08_glycolysis_gradient.R`.
- Figure 9B: spatial GlycoHigh versus GlycoLow direction consistency.
- Supplementary Figure S20: RCTD-estimated spatial cell-type composition.
- Supplementary Figure S21: metabolic specificity / partial Spearman correlation.
- Supplementary Figure S22: NicheNet ligand activity.
- Supplementary Figure S23: patient-level ligand mean-expression robustness.
- Supplementary Figure S24: TCGA-LIHC and GSE14520 3-year OS calibration.
- Supplementary Figure S25: GCK sensitivity analysis.

Tentative or scope-pending figures are not included in the final restoration scripts until their data scope is finalized.

## Data

- scRNA-seq: GSE149614 (primary), GSE125449 (validation)
- Bulk RNA-seq: TCGA-LIHC (UCSC Xena)
- Microarray: GSE14520
- Spatial transcriptomics: GSE238264
