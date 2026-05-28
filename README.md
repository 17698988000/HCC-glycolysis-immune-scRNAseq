# HCC glycolysis-immune scRNA-seq repository

This repository contains analysis and reproducibility materials for the current JTM revision of the HCC GlycoHigh malignant-hepatocyte metabolic-immune state manuscript.

**Current manuscript alignment:** Step-8-final claim-control version (`JTM_revision_working_version_v1_2026-05-27_step8_claim_control_final.docx`).

## Interpretation boundary

This repository supports a **discovery-stage / hypothesis-generating** transcriptomic study. It should not be read as a clinical-ready prognostic assay, a treatment-selection model, or a mechanistic validation package.

- The four-gene score is a `TPI1/ENO1/LDHA/SLC2A1` **tissue-level readout**, not an implemented clinical assay.
- `SPP1` and `MIF` analyses support candidate / inferred / complementary communication or transcriptional-response layers only.
- `GSE235863` is used only for exploratory anti-PD-1 plus lenvatinib association analysis (`n = 15`; non-responders `n = 4`). The active repository does not include ROC/PR, DCA, complex ML, or treatment-selection analysis for this cohort.
- inferCNV outputs are used as malignancy/CNV support only, not as mechanistic validation.

## Locked analysis conventions

- Primary object: `seurat_final.rds`
- Patient column: `patient`
- Cell-type column: `cell_type`
- Site column: `site`
- Glycolysis activity column: `Glycolysis_AUC`
- Tumor-derived hepatocytes: `site == "Tumor" & cell_type == "Hepatocyte"`
- Tumor-derived hepatocytes: `n = 15,391`
- GlycoHigh cells: `n = 7,695`
- GlycoLow cells: `n = 7,696`

The locked TCGA-derived four-gene tissue-level readout is:

```text
score = 0.3041908 * TPI1 +
        0.9639654 * ENO1 +
        1.3404374 * LDHA +
        0.2424239 * SLC2A1
```

The score was derived from TCGA-LIHC primary tumor survival-matched samples (`n = 365`). The `lambda.min` model retained four genes; the conservative `lambda.1se` model was null and is treated as an important limitation.

## Manuscript-to-repository map

Use `MANUSCRIPT_TO_REPOSITORY_MAP.tsv` as the authoritative map between the current manuscript, supplementary materials, scripts, and interpretation boundaries.

## Recommended run order

Use `RUN_ORDER.md`. Some large public-data downloads, Seurat object construction, inferCNV objects, and locked figure outputs require local source objects or previously generated intermediate files. Scripts are intended for transparent reproducibility and source/QC generation, not for changing the locked manuscript figures unless a new QC pass is documented.

## Active script status

| File | Current role |
|---|---|
| `01_QC_clustering.R` | Upstream single-cell QC, integration, clustering, and annotation workflow. |
| `02_glycolysis_scoring.R` | Locked GlycoHigh/GlycoLow source/QC generation using current metadata and `Glycolysis_AUC`. |
| `03_cellchat_analysis.R` | CellChat source/QC workflow; inference-based communication support only. |
| `06_inferCNV_malignant.R` | Marker-score source/QC and optional inferCNV workflow for malignancy/CNV support. |
| `07_LASSO_risk_score.R` | TCGA-derived four-gene LASSO/Cox tissue-level readout construction. |
| `08_glycolysis_gradient.R` | Auxiliary glycolysis-gradient source/QC analysis in tumor-derived hepatocytes. |
| `09_TF_activity.R` | Auxiliary TF-activity source/QC analysis; not a mechanistic validation assay. |
| `10_GSE14520_validation.R` | External score-level survival evaluation of the locked four-gene readout in GSE14520. |
| `11_spatial_transcriptomics.R` | GSE238264 Visium spatial analysis and RCTD-related source/QC outputs. |
| `12_GSE235863_exploratory_association.R` | Current GSE235863 descriptive, hypothesis-generating association script. No ROC/PR/DCA/ML. |
| `13_revision_analyses.R` | Status/index script only; not a source of additional biological or clinical claims. |
| `14_OXPHOS_metabolic_specificity.R` | OXPHOS/metabolic-specificity source/QC support. |
| `15_NicheNet_analysis.R` | NicheNet locked-status / reproduction-guard script. Inference-based only. |
| `16_patient_level_ligand_effects.R` | Patient-level ligand mean-expression robustness source/QC analysis. |
| `17_partial_correlation_metabolic_specificity.R` | Partial Spearman metabolic-specificity analysis for glycolysis/OXPHOS context. |
| `revision_figure_restore/` | Selected vector-output restoration scripts retained for current supplementary/model-diagnostic support. |
| `archive_not_used/` | Scripts retained for transparency but not used in the current Step-8-final manuscript. |

## Archived / not used in the current manuscript

The following were moved out of the active root directory because they conflict with the current manuscript scope or were removed during earlier revision steps:

- `05_GSE125449_validation.R`: GSE125449 was removed from the current manuscript dataset census and should not appear as an active dataset.
- `12_drug_repurposing.R`: legacy GSE235863 performance-oriented script. Replaced by `12_GSE235863_exploratory_association.R`.
- ENO1-centered TCGA and figure-restoration scripts: retained only for historical traceability. The current manuscript is not an ENO1 single-gene story.

See `archive_not_used/README.md` for details.

## Data resources

Public datasets used in the current manuscript:

- scRNA-seq discovery: `GSE149614`
- Spatial transcriptomics: `GSE238264`
- Bulk tissue-level discovery / clinical association: `TCGA-LIHC`
- External score-level survival evaluation: `GSE14520`
- Exploratory anti-PD-1 plus lenvatinib cohort: `GSE235863`

`GSE125449` is archived/not used in the current manuscript version.

## Reproducibility notes

- This repository does not upload protected or restricted raw data.
- Some scripts require locally downloaded public data or intermediate R objects.
- Use source/QC outputs and the manifest files generated by each script to document local reruns.
- Do not overwrite locked manuscript figures without updating the corresponding QC outputs, figure legends, and manuscript-to-repository map.
- README text is not a substitute for running the script-level QC checks.

## Code availability note

The manuscript currently contains a placeholder `[GitHub URL]` in the Code availability section. Replace it with the final public repository URL before submission.
