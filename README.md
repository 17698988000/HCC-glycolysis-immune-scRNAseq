# HCC-glycolysis-immune-scRNAseq

Analysis code for the ENO1/glycolysis-immune evasion study in hepatocellular carcinoma (HCC).

This repository contains the analysis scripts used for single-cell, spatial, bulk-transcriptomic, and validation analyses. Revision-stage scripts have been converted to source/QC-first or locked-status workflows to prevent obsolete local-path code, outdated cell labels, or changed object definitions from regenerating manuscript-inconsistent outputs.

## Locked analysis conventions

The current single-cell analysis convention is:

- Primary object: `seurat_final.rds`
- Patient column: `patient`
- Cell-type column: `cell_type`
- Site column: `site`
- Glycolysis activity column: `Glycolysis_AUC`
- Tumor-derived hepatocytes: `site == "Tumor" & cell_type == "Hepatocyte"`
- Number of tumor-derived hepatocytes: `n = 15,391`
- Glycolysis median cutoff: `0.2203849`
- GlycoHigh definition: `Glycolysis_AUC > median_cut`
- GlycoLow definition: remaining tumor-derived hepatocytes
- GlycoLow cells: `7,696`
- GlycoHigh cells: `7,695`

The locked TCGA-derived four-gene risk score is:

```text
score = 0.3041908 * TPI1 +
        0.9639654 * ENO1 +
        1.3404374 * LDHA +
        0.2424239 * SLC2A1
```

The four-gene score was derived from TCGA-LIHC primary tumor survival-matched samples with `n = 365`. The `lambda.min` model retained four genes, whereas the `lambda.1se` model was null.

## Current script status

| Script | Current role / status |
|---|---|
| `01_QC_clustering.R` | Primary scRNA-seq quality control, integration, clustering, and annotation workflow. Leave unchanged unless the upstream object-building workflow is intentionally being re-derived. |
| `02_glycolysis_scoring.R` | Locked current-object GlycoHigh/GlycoLow assignment and QC script. It reads `seurat_final.rds`, uses existing `Glycolysis_AUC`, validates `patient`, `cell_type`, `site`, and `Glycolysis_AUC`, enforces tumor-derived hepatocyte `n = 15,391`, median cutoff `0.2203849`, GlycoLow `n = 7,696`, and GlycoHigh `n = 7,695`, then writes source assignment/QC outputs. It does not recalculate AUCell or overwrite the primary `seurat_final.rds`. |
| `03_cellchat_analysis.R` | Updated current-object CellChat source/QC workflow. Legacy `D:/scRNA_project` paths, `tumor_hepatocytes.rds` dependency, and obsolete `Hepatocytes` label have been removed. Round 1 uses original `cell_type`; Round 2 uses GlycoHigh/GlycoLow tumor hepatocytes plus remaining cell types; Round 3 is optional only when an explicit subtype column already exists. |
| `04_TCGA_validation.R` | TCGA-LIHC ENO1 survival, immune infiltration, and TIDE-related validation analyses. |
| `05_GSE125449_validation.R` | Updated source/QC-first GSE125449 validation workflow. Legacy local paths and fixed hepatocyte-cluster assumptions have been removed. The script loads a prebuilt object or Set1/Set2 matrices, identifies hepatocyte-like cells using available annotation or marker scoring, scores glycolysis, and writes validation source/QC outputs. |
| `06_inferCNV_malignant.R` | Updated current-object malignant hepatocyte marker-score/QC workflow with optional inferCNV. Legacy local paths, automatic gene-order download, and obsolete `Hepatocytes` label requirement have been removed. Optional inferCNV runs only if the package and a local `hg38_gencode_v27.txt` gene-order file are available. |
| `07_LASSO_risk_score.R` | Locked TCGA-derived four-gene Cox/LASSO risk-score model using `TPI1`, `ENO1`, `LDHA`, and `SLC2A1`. |
| `08_glycolysis_gradient.R` | Supplementary Figure S15 glycolysis-gradient analysis using 10 equal-width AUCell-score bins and mean log-normalized expression of `ENO1`, `LDHA`, `SPP1`, `MIF`, and `PTGES` in tumor-derived hepatocytes. |
| `09_TF_activity.R` | Supplementary Figure S16 DoRothEA/VIPER transcription-factor activity analysis in 15,391 tumor-derived hepatocytes, correlated with `Glycolysis_AUC`. |
| `10_GSE14520_validation.R` | Independent GSE14520 validation of the locked four-gene score. Tumor samples with complete survival information: `n = 221`; multivariable Cox evaluable samples: `n = 217`; HR per SD approximately `1.32`. |
| `11_spatial_transcriptomics.R` | GSE238264 Visium spatial analysis for HCC1R, HCC2R, HCC3R, and HCC4R using SCTransform and within-sample median Glycolysis1 module-score grouping. |
| `12_drug_repurposing.R` | Replaced legacy drug-repurposing enrichment with GSE235863 anti-PD-1 plus lenvatinib exploratory non-response association analysis using the locked four-gene score. This is exploratory and not a formal response-prediction model. |
| `13_revision_analyses.R` | Legacy/status note script. It no longer regenerates mixed obsolete revision figures and redirects analyses to dedicated current scripts. |
| `14_OXPHOS_metabolic_specificity.R` | OXPHOS AUCell score and S21 input source/QC generation for tumor-derived hepatocytes. This script does not generate S22 and does not generate the final S21 figure. |
| `15_NicheNet_analysis.R` | Supplementary Figure S22 locked-status / reproduction-guard script. It intentionally does not rerun NicheNet and does not overwrite the manuscript-matched S22 figure. |
| `16_patient_level_ligand_effects.R` | Supplementary Figure S23 patient-level GlycoHigh-minus-GlycoLow mean-expression robustness analysis for `SPP1`, `MIF`, and `PTGES`. |
| `17_partial_correlation_metabolic_specificity.R` | Supplementary Figure S21 partial Spearman metabolic-specificity analysis for glycolysis versus OXPHOS, using `SPP1` and `MIF`. |
| `restore_FigS12_ENO1_glycolysis_per_patient.R` | Final restoration script for Supplementary Figure S12 per-patient ENO1-glycolysis correlation. |
| `revision_figure_restore/` | Final submission-ready clean vector figure restoration scripts after figure-level QC, renumbering, and vector-output restoration. |

## Early-script repair status

The early single-cell scripts have the following current maintenance status:

| Script | Repair status |
|---|---|
| `02_glycolysis_scoring.R` | Current and retained. It is a locked assignment/QC script based on existing `Glycolysis_AUC`. Do not replace it with a re-scoring script unless a full AUCell recomputation is intentionally required. |
| `03_cellchat_analysis.R` | Repaired. Uses `seurat_final.rds` and locked metadata columns; exports CellChat source/QC outputs rather than direct manuscript figures. |
| `05_GSE125449_validation.R` | Repaired. Removes local paths and fixed hepatocyte cluster IDs; uses annotation or marker scoring to identify hepatocyte-like cells. |
| `06_inferCNV_malignant.R` | Repaired. Removes local paths and obsolete labels; performs marker-score QC and optional inferCNV only when required local dependencies exist. |
| `01_QC_clustering.R` | Not revised in this round. Leave as the upstream object-building script unless the full object-generation workflow is being re-audited. |

## Source/QC-first figure policy

For revised scripts, the intended workflow is:

1. Write source-data tables.
2. Write QC tables.
3. Generate PDF/PNG figure outputs only after mandatory QC checks pass.
4. Avoid in-figure titles such as “Supplementary Figure Sxx” in final figure panels unless explicitly required.
5. Do not regenerate locked manuscript figures from obsolete local paths or changed object definitions.

## Locked and corrected figure mapping

| Figure | Current source / status |
|---|---|
| Supplementary Figure S12 | Per-patient ENO1-glycolysis correlation in tumor-derived hepatocytes; restored by `restore_FigS12_ENO1_glycolysis_per_patient.R`. |
| Supplementary Figure S15 | Glycolysis-gradient expression of `ENO1`, `LDHA`, `SPP1`, `MIF`, and `PTGES` in tumor-derived hepatocytes; reproduced by `08_glycolysis_gradient.R`. |
| Supplementary Figure S16 | DoRothEA/VIPER transcription-factor activity correlated with `Glycolysis_AUC`; reproduced by `09_TF_activity.R`. |
| Supplementary Figure S18 | Spatial transcriptomics output from `11_spatial_transcriptomics.R`. |
| Supplementary Figure S20 | RCTD-estimated spatial cell-type composition from `11_spatial_transcriptomics.R`. |
| Supplementary Figure S21 | Glycolysis-versus-OXPHOS partial Spearman metabolic-specificity analysis; source/QC from `14_OXPHOS_metabolic_specificity.R`, final plot from `17_partial_correlation_metabolic_specificity.R`. |
| Supplementary Figure S22 | NicheNet ligand activity analysis; locked manuscript-matched status protected by `15_NicheNet_analysis.R`. |
| Supplementary Figure S23 | Patient-level ligand mean-expression robustness for `SPP1`, `MIF`, and `PTGES`; reproduced by `16_patient_level_ligand_effects.R`. |
| Supplementary Figure S24 | TCGA-LIHC and GSE14520 3-year OS calibration. |
| Supplementary Figure S25 | GCK sensitivity analysis. |

Tentative or scope-pending figures are not included in final restoration scripts until their data scope is finalized.

## Locked supplementary results

### Supplementary Figure S12

S12 is locked to the current-object result:

- Tumor-derived hepatocytes: `15,391`
- HCC patients: `8`
- Global Pearson correlation: `R = 0.57`
- Per-patient Spearman rho range: `0.363` in HCC04 to `0.633` in HCC10
- Median per-patient rho: `0.498`
- Raw p-value range: `4.43e-248` to `5.06e-17`
- All BH-adjusted p-values `< 0.05`

### Supplementary Figure S15

S15 is locked as a glycolysis-gradient analysis:

- Tumor-derived hepatocytes only
- 10 equal-width `Glycolysis_AUC` bins
- Mean log-normalized expression of `ENO1`, `LDHA`, `SPP1`, `MIF`, and `PTGES`
- Reproduced by `08_glycolysis_gradient.R`

### Supplementary Figure S21

S21 is the metabolic-specificity / partial Spearman analysis:

- `SPP1` glycolysis partial rho approximately `0.227`
- `SPP1` OXPHOS partial rho approximately `0.125`
- `MIF` glycolysis partial rho approximately `0.105`
- `MIF` OXPHOS partial rho approximately `0.242`
- Glycolysis-OXPHOS rho approximately `0.066`

### Supplementary Figure S22

S22 is locked to the manuscript-matched NicheNet ligand-activity result:

- Sender cells: GlycoHigh hepatocytes, `n = 5,589`
- Receiver cells: tumor-derived T/NK plus myeloid cells, `n = 11,383`
- Expressed ligands: `325`
- Ligand activity panel: top 30 ligands
- `MIF` rank: `37`
- `MIF` AUPR: `0.091`
- `SPP1` rank: `250`
- `SPP1` AUPR: `0.026`

`15_NicheNet_analysis.R` intentionally does not rerun NicheNet because regenerating S22 from the current `seurat_final.rds` median split would use GlycoHigh `n = 7,695`, which may not reproduce the manuscript-matched S22 result. If a full S22 computational rerun is required, the original S22-specific NicheNet input/output objects that produced sender `n = 5,589` should be restored first.

### Supplementary Figure S23

S23 is locked to the patient-level ligand mean-expression robustness result:

- `MIF` mean effect High-Low: `0.460`
- `MIF` median effect High-Low: `0.355`
- `MIF` paired Wilcoxon p-value: `0.0143`
- `MIF` paired Wilcoxon FDR: `0.0428`
- `SPP1` paired Wilcoxon FDR: `0.441`
- `PTGES` paired Wilcoxon FDR: `0.0888`
- Locked final files: `FigS23_final_no_title.pdf` and `FigS23_final_no_title.png`

## Spatial transcriptomics notes

`11_spatial_transcriptomics.R` uses GSE238264 Visium samples:

- HCC1R
- HCC2R
- HCC3R
- HCC4R

Expected within-sample median group counts are:

| Sample | GlycoHigh | GlycoLow |
|---|---:|---:|
| HCC1R | 1,503 | 1,503 |
| HCC2R | 1,383 | 1,383 |
| HCC3R | 1,085 | 1,085 |
| HCC4R | 1,501 | 1,501 |

Spatial conclusions should be described as spot-level tissue co-enrichment. They should not be overstated as same-cell co-expression, tumor-cell-specific ligand production, or direct ligand-receptor contact.

## GSE125449 validation notes

`05_GSE125449_validation.R` is now a source/QC-first validation workflow. It should be interpreted as a validation analysis rather than a manuscript-locked figure restoration script.

Key safeguards:

- No `D:/scRNA_project` local path is required.
- Fixed hepatocyte cluster IDs are not used.
- Existing annotation is preferred when available.
- Marker scoring is used only when a suitable annotation column is absent.
- GlycoHigh/GlycoLow is defined by median split within hepatocyte-like cells in the validation dataset.

## inferCNV / malignant hepatocyte notes

`06_inferCNV_malignant.R` is now a current-object source/QC workflow with optional inferCNV.

Key safeguards:

- Tumor-derived hepatocytes are selected using `site == "Tumor" & cell_type == "Hepatocyte"`.
- The locked expected tumor-derived hepatocyte count is `15,391`.
- The script writes marker-score source/QC outputs even when inferCNV is not available.
- inferCNV runs only if the `infercnv` package and a local `hg38_gencode_v27.txt` gene-order file are available.
- The script does not automatically download gene-order files.

## GSE235863 exploratory response analysis

`12_drug_repurposing.R` currently represents an exploratory anti-PD-1 plus lenvatinib non-response association analysis, not a drug-repurposing prioritization analysis.

Locked exploratory values:

- Total samples: `n = 15`
- Responders: `11`
- Non-responders: `4`
- Positive class: non-responder
- Four-gene score AUC approximately `0.932`
- 95% CI approximately `0.773–1.000`
- Median-in-high grouping places all 4 non-responders in the High group
- Fisher p-value approximately `0.077`

## Data resources

- scRNA-seq: GSE149614 primary cohort, GSE125449 validation cohort
- Bulk RNA-seq: TCGA-LIHC
- Microarray validation: GSE14520
- Spatial transcriptomics: GSE238264
- Exploratory treatment-response association: GSE235863

## Notes for future maintenance

- Update `README.md` only after script-level changes are committed.
- Do not use README text as evidence that a script has already been committed; verify the actual file on `main`.
- For locked figures, do not replace manuscript-matched results unless the exact source objects and QC targets are available.
- Legacy or pending-review scripts should not be used to overwrite locked figure outputs without a new, documented QC pass.
- Do not replace the current `02_glycolysis_scoring.R` with a re-scoring/recalculation script unless the goal is an explicit full AUCell recomputation and new downstream QC.
