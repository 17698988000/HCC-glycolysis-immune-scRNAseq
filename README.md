# HCC-glycolysis-immune-scRNAseq

Analysis code for the manuscript:

**Integrated single-cell, spatial, and clinical transcriptomics delineate a glycolysis-high metabolic-immune refractory subtype in hepatocellular carcinoma**

This repository contains scripts for the single-cell, spatial, and clinical transcriptomic analyses used to define a glycolysis-high metabolic-immune refractory subtype in hepatocellular carcinoma.

## Current finalized data convention

The finalized single-cell analysis uses the current integrated Seurat object:

- Main input object: `seurat_final.rds`
- Patient column: `patient`
- Cell-type column: `cell_type`
- Tissue/site column: `site`
- Glycolysis AUCell column: `Glycolysis_AUC`

Tumor-derived hepatocytes are defined as:

    site == "Tumor" & cell_type == "Hepatocyte"

Locked tumor-derived hepatocyte count:

    n = 15,391

GlycoHigh and GlycoLow are defined by a median split of `Glycolysis_AUC` within tumor-derived hepatocytes only:

    median_cut <- median(Glycolysis_AUC)
    High <- Glycolysis_AUC > median_cut
    Low  <- Glycolysis_AUC <= median_cut

Locked values:

    Median Glycolysis_AUC cutoff = 0.2203849
    GlycoLow  = 7,696
    GlycoHigh = 7,695

Deprecated conventions that should not be used in final manuscript-aligned scripts:

- `glycolysis_score` instead of `Glycolysis_AUC`
- `tissue` instead of `site`
- `Hepatocytes` instead of `Hepatocyte`
- Stratifying all cells instead of tumor-derived hepatocytes
- Using `>= median` for High; final High is strictly `> median`

## Scripts

| Script | Description |
|---|---|
| `01_QC_clustering.R` | Quality control, Harmony batch correction, and cell-type annotation. |
| `02_glycolysis_scoring.R` | Finalized tumor-derived hepatocyte GlycoHigh/GlycoLow assignment using `seurat_final.rds`, `Glycolysis_AUC`, `site == "Tumor"`, `cell_type == "Hepatocyte"`, and locked QC checks for n = 15,391, GlycoLow = 7,696, GlycoHigh = 7,695. |
| `03_cellchat_analysis.R` | CellChat communication analysis comparing GlycoHigh and GlycoLow tumor-derived hepatocytes with immune target populations. |
| `04_TCGA_validation.R` | TCGA-LIHC ENO1 survival, immune infiltration, and TIDE-signature-derived immune-context analyses. |
| `05_GSE125449_validation.R` | Cross-dataset evaluation in GSE125449 for ENO1-glycolysis association and ligand expression. |
| `06_inferCNV_malignant.R` | inferCNV-based malignant hepatocyte support analysis; hepatocyte labels should be aligned with the finalized `Hepatocyte` convention. |
| `07_LASSO_risk_score.R` | TCGA-LIHC four-gene glycolysis score construction using TPI1, ENO1, LDHA, and SLC2A1 with locked coefficients; lambda.min retains four genes and lambda.1se is null. |
| `08_glycolysis_gradient.R` | Supplementary Figure S15 glycolysis-gradient analysis using 10 equal-width AUCell-score bins and mean log-normalized expression of ENO1, LDHA, SPP1, MIF, and PTGES in tumor-derived hepatocytes. |
| `09_TF_activity.R` | DoRothEA/VIPER transcription factor activity inference in tumor-derived hepatocytes using `Glycolysis_AUC`. |
| `10_GSE14520_validation.R` | External validation of the finalized four-gene glycolysis score in GSE14520; survival cohort n = 221 and multivariable Cox n = 217. |
| `11_spatial_transcriptomics.R` | Visium spatial transcriptomics analysis of glycolysis and immunosuppressive ligand features. |
| `12_drug_repurposing.R` | Legacy exploratory drug-repurposing script; not part of the finalized main manuscript line unless replaced by the GSE235863 immunotherapy non-response analysis. |
| `13_revision_analyses.R` | Revision-stage supplementary analyses and figure restoration support. |
| `14_OXPHOS_metabolic_specificity.R` | OXPHOS AUCell scoring and metabolic specificity analysis; final output corresponds to Supplementary Figure S21. |
| `15_NicheNet_analysis.R` | NicheNet ligand activity analysis using GlycoHigh tumor-derived hepatocytes as sender cells and tumor-derived T/NK plus myeloid cells as receiver cells. |
| `16_patient_level_ligand_effects.R` | Patient-level ligand mean-expression robustness analysis for SPP1, MIF, and PTGES; final MIF paired Wilcoxon FDR = 0.0428. |
| `17_partial_correlation_metabolic_specificity.R` | Supplementary Figure S21 partial Spearman correlation plot for glycolysis versus OXPHOS specificity in tumor-derived hepatocytes. |

## Corrected final figure mapping

- Supplementary Figure S12: per-patient ENO1-glycolysis correlation in tumor-derived hepatocytes, reproduced using the current-object `seurat_final.rds` convention.
- Supplementary Figure S15: glycolysis-gradient expression of ENO1, LDHA, SPP1, MIF, and PTGES in tumor-derived hepatocytes, reproduced by `08_glycolysis_gradient.R`.
- Supplementary Figure S20: RCTD-estimated spatial cell-type composition.
- Supplementary Figure S21: metabolic specificity / partial Spearman correlation.
- Supplementary Figure S22: NicheNet ligand activity.
- Supplementary Figure S23: patient-level ligand mean-expression robustness.
- Supplementary Figure S24: TCGA-LIHC and GSE14520 3-year OS calibration.
- Supplementary Figure S25: GCK inclusion/exclusion sensitivity analysis.

## Locked finalized supplementary figure results

### Supplementary Figure S12

- Input object: `seurat_final.rds`
- Tumor-derived hepatocytes: 15,391
- Patients: 8 HCC patients
- Global Pearson R = 0.57
- Per-patient Spearman rho range = 0.363 in HCC04 to 0.633 in HCC10
- Median rho = 0.498
- Raw p-value range = 4.43e-248 to 5.06e-17
- All BH-adjusted p-values < 0.05
- Final files: `FigS12_final_two_color_no_stars.pdf`, `FigS12_final_two_color_no_stars.png`

### Supplementary Figure S15

- Input object: `seurat_final.rds`
- Tumor-derived hepatocytes: 15,391
- Glycolysis column: `Glycolysis_AUC`
- Binning method: 10 equal-width bins spanning the AUCell glycolysis score range
- Binning code: `cut(Glycolysis_AUC, breaks = 10, include.lowest = TRUE)`
- Expression summary: mean log-normalized RNA expression
- Genes: ENO1, LDHA, SPP1, MIF, PTGES
- Final output files:
  - `FigS15_source_data.csv`
  - `FigS15_QC_check.csv`
  - `FigS15_final_no_title.pdf`
  - `FigS15_final_no_title.png`

### Supplementary Figure S23

- MIF mean effect High-Low = 0.460
- MIF median effect High-Low = 0.355
- MIF paired Wilcoxon p = 0.0143
- MIF paired Wilcoxon FDR = 0.0428
- SPP1 paired Wilcoxon FDR = 0.441
- PTGES paired Wilcoxon FDR = 0.0888
- Final files: `FigS23_final_no_title.pdf`, `FigS23_final_no_title.png`

## Four-gene glycolysis score

The finalized tissue-level glycolysis score uses four genes:

    TPI1, ENO1, LDHA, SLC2A1

Score formula:

    Four-gene glycolysis score =
    0.3041908 * TPI1 +
    0.9639654 * ENO1 +
    1.3404374 * LDHA +
    0.2424239 * SLC2A1

The old nine-gene score is deprecated and should not be used in final manuscript-aligned scripts.

TCGA-LIHC finalized LASSO convention:

- Primary tumor survival-matched cohort: n = 365
- Candidate features: 22 curated glycolysis genes
- `lambda.min`: retained TPI1, ENO1, LDHA, and SLC2A1
- `lambda.1se`: null model

GSE14520 finalized validation convention:

- Tumor samples with survival data: n = 221
- Multivariable Cox complete-case cohort: n = 217
- Adjusted covariates:
  - AFP
  - cirrhosis
  - main tumor size
  - multinodular disease

## Data

- scRNA-seq: GSE149614 primary dataset, GSE125449 cross-dataset evaluation
- Bulk RNA-seq: TCGA-LIHC from UCSC Xena
- Microarray: GSE14520
- Spatial transcriptomics: GSE238264
- Exploratory immunotherapy cohort: GSE235863

## Repository status note

This repository is being synchronized to the finalized manuscript. Scripts that still contain legacy assumptions should be updated before being treated as manuscript-reproducible.

Priority corrections include:

1. `02_glycolysis_scoring.R`: finalized `Glycolysis_AUC` and tumor-derived hepatocyte split.
2. `07_LASSO_risk_score.R`: finalized four-gene TCGA score.
3. `10_GSE14520_validation.R`: finalized four-gene GSE14520 validation.
4. `12_drug_repurposing.R`: legacy status or replacement by GSE235863 non-response analysis.
5. `14_OXPHOS_metabolic_specificity.R`: final Supplementary Figure S21 mapping.
6. `15_NicheNet_analysis.R`: current-object and site/cell-type convention alignment.
7. `17_partial_correlation_metabolic_specificity.R`: tumor-derived hepatocyte restriction, n = 15,391.
