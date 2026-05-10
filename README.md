
# HCC-glycolysis-immune-scRNAseq

Analysis code for ENO1/glycolysis-immune evasion study in HCC

## Scripts

| Script | Description |
|--------|-------------|
| 01_QC_clustering.R | Quality control, Harmony batch correction, cell type annotation (Sections 2.1–2.3) |
| 02_glycolysis_scoring.R | AUCell glycolysis scoring, GlycoHigh/GlycoLow stratification (Section 2.4) |
| 03_cellchat_analysis.R | CellChat cell-cell communication analysis (Section 2.6) |
| 04_TCGA_validation.R | TCGA-LIHC ENO1 survival analysis, immune infiltration, TIDE scores (Sections 2.7–2.8) |
| 05_GSE125449_validation.R | Cross-dataset validation in GSE125449 (Section 3.8) |
| 06_inferCNV_malignant.R | inferCNV malignant hepatocyte identification (Section 2.5) |
| 07_LASSO_risk_score.R | LASSO-penalized Cox regression glycolysis risk score (Section 2.10) |
| 08_glycolysis_gradient.R | Glycolysis activity gradient analysis (Section 2.11) |
| 09_TF_activity.R | DoRothEA/VIPER transcription factor activity inference (Section 2.12) |
| 10_GSE14520_validation.R | Independent validation of risk score in GSE14520 (Section 2.13) |
| 11_spatial_transcriptomics.R | Visium spatial transcriptomics + RCTD deconvolution (Section 2.14) |
| 12_drug_repurposing.R | Drug repurposing enrichment analysis (Section 2.17) |
| 13_revision_analyses.R | Revision-stage supplementary analyses (Figures S9–S13, S20) |
| 14_OXPHOS_metabolic_specificity.R | OXPHOS AUCell scoring and partial Spearman correlations (Section 2.15, Figure S22) |
| 15_NicheNet_analysis.R | NicheNet ligand activity analysis (Section 2.16, Figure S23) |

## Data

- scRNA-seq: GSE149614 (primary), GSE125449 (validation)
- Bulk RNA-seq: TCGA-LIHC (UCSC Xena)
- Microarray: GSE14520
- Spatial transcriptomics: GSE238264
