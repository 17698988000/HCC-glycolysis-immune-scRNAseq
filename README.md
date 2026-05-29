# HCC GlycoHigh malignant-hepatocyte state repository

This repository contains analysis scripts and reproducibility documentation for the revised Journal of Translational Medicine manuscript on a GlycoHigh malignant-hepatocyte metabolic-immune state and a compact four-gene tissue-level readout in hepatocellular carcinoma.

**Current manuscript alignment:** current revised manuscript package  
`JTM_revision_working_version_v1_2026-05-27_step11_crosscheck_minorfix_v6.docx`

## Study scope and interpretation boundaries

This repository supports a discovery-stage and hypothesis-generating transcriptomic study. It should not be interpreted as a clinical-ready prognostic assay, treatment-selection model, or mechanistic validation package.

Key boundaries:

* The TPI1/ENO1/LDHA/SLC2A1 score is a tissue-level readout of GlycoHigh biology, not an implemented clinical assay.
* ENO1 is retained as an interpretable component of the four-gene readout, not as a stand-alone mechanistic driver or single-gene main theme.
* SPP1 and MIF analyses support candidate, inferred, complementary communication or transcriptional-response layers only.
* CellChat, NicheNet, expression-level ligand-receptor compatibility, and LIANA analyses are inference-based transcriptomic support, not functional validation of ligand secretion, receptor activation, immune suppression, or treatment response.
* GSE235863 is used only for exploratory anti-PD-1 plus lenvatinib association analysis (n = 15; non-responders, n = 4). The active repository does not include ROC, PR, DCA, complex machine learning, or treatment-selection analysis for this cohort.
* inferCNV outputs are used as malignancy/CNV support only, not as mechanistic validation.
* GSE125449 is archived/not used in the current manuscript package.

## Locked analysis conventions

The locked single-cell analysis conventions are:

* Primary object: `seurat_final.rds`
* Patient column: `patient`
* Cell-type column: `cell_type`
* Site column: `site`
* Glycolysis activity column: `Glycolysis_AUC`
* Tumor-derived hepatocytes: `site == "Tumor" & cell_type == "Hepatocyte"`
* Tumor-derived hepatocytes: n = 15,391
* GlycoHigh cells: n = 7,695
* GlycoLow cells: n = 7,696

The locked TCGA-derived four-gene tissue-level readout is:

```text
score = 0.3041908 * TPI1 +
        0.9639654 * ENO1 +
        1.3404374 * LDHA +
        0.2424239 * SLC2A1
```

The score was derived from TCGA-LIHC primary tumor survival-matched samples (n = 365). The lambda.min model retained four genes; the conservative lambda.1se model was null and is treated as an important limitation.

## Data resources

Public datasets used in the current manuscript package:

| Dataset   | Role in manuscript                                   | Boundary                                                                                    |
| --------- | ---------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| GSE149614 | Single-cell discovery layer                          | Paired tumor/non-tumor HCC single-cell analysis; eight paired patients used after filtering |
| GSE238264 | Spatial transcriptomic support                       | Four Visium HCC sections; spot-level co-enrichment only                                     |
| TCGA-LIHC | Bulk tissue-level discovery and clinical association | Retrospective score-level association; not clinical deployment                              |
| GSE14520  | External score-level survival evaluation             | External survival association using fixed TCGA-derived coefficients                         |
| GSE235863 | Exploratory anti-PD-1 plus lenvatinib association    | Small, class-imbalanced cohort; hypothesis-generating only                                  |
| GSE125449 | Archived/not used                                    | Not part of the active current manuscript package                                           |

## Active scripts

| File                                             | Current role                                                                                             |
| ------------------------------------------------ | -------------------------------------------------------------------------------------------------------- |
| `01_QC_clustering.R`                             | Upstream single-cell QC, integration, clustering, and annotation workflow                                |
| `02_glycolysis_scoring.R`                        | GlycoHigh/GlycoLow source and QC generation using the locked glycolysis activity score                   |
| `03_cellchat_analysis.R`                         | CellChat source/QC workflow; inference-based communication support only                                  |
| `04_TCGA_validation.R`                           | TCGA-LIHC survival and clinical association support                                                      |
| `06_inferCNV_malignant.R`                        | Marker-score and inferCNV/CNV support for malignant-feature assessment                                   |
| `07_LASSO_risk_score.R`                          | TCGA-derived four-gene LASSO/Cox tissue-level readout construction                                       |
| `08_glycolysis_gradient.R`                       | Auxiliary glycolysis-gradient source/QC analysis in tumor-derived hepatocytes                            |
| `09_TF_activity.R`                               | Auxiliary TF-activity source/QC analysis; not a mechanistic validation assay                             |
| `10_GSE14520_validation.R`                       | External score-level survival evaluation of the locked four-gene readout in GSE14520                     |
| `11_spatial_transcriptomics.R`                   | GSE238264 Visium spatial analysis and RCTD-related source/QC outputs                                     |
| `12_GSE235863_exploratory_association.R`         | GSE235863 descriptive, hypothesis-generating association script; no ROC/PR/DCA/ML                        |
| `13_revision_analyses.R`                         | Status/index script only; not a source of additional biological or clinical claims                       |
| `14_OXPHOS_metabolic_specificity.R`              | OXPHOS/metabolic-specificity source/QC support                                                           |
| `15_NicheNet_analysis.R`                         | NicheNet reproduction-guard and sensitivity script; inference-based only                                 |
| `16_patient_level_ligand_effects.R`              | Patient-level ligand mean-expression robustness source/QC analysis                                       |
| `17_partial_correlation_metabolic_specificity.R` | Partial Spearman metabolic-specificity analysis for glycolysis/OXPHOS context                            |
| `18_LIANA_directional_concordance.R`             | LIANA directional-concordance check for predefined MIF/SPP1 axes; candidate inference-based support only |

## Archived and historical scripts

The `archive_not_used/` folder contains scripts retained for transparency but not used in the current revised manuscript package.

The following analyses are not part of the active manuscript workflow:

* GSE125449 validation scripts
* Legacy GSE235863 ROC/AUC, drug-repurposing, or treatment-selection scripts
* ENO1-centered restoration scripts
* Historical figure-restoration scripts that are not authoritative for the current Figure 1-7 workflow

Archived scripts should not be cited as evidence for the submitted manuscript.

## Manuscript-to-repository map

Use `MANUSCRIPT_TO_REPOSITORY_MAP.tsv` as the authoritative map between the current manuscript, supplementary materials, scripts, and interpretation boundaries.

Use `RUN_ORDER.md` for the recommended execution order. Some large public-data downloads, Seurat object construction, inferCNV objects, RCTD objects, and locked figure outputs require local source objects or previously generated intermediate files. Scripts are intended for transparent reproducibility and source/QC generation, not for changing locked manuscript figures unless a new QC pass is documented.

## Code availability

All R analysis scripts used in this study are available at:

```text
https://github.com/17698988000/HCC-glycolysis-immune-scRNAseq
```

The repository supports the current discovery-stage transcriptomic manuscript package and should not be interpreted as a clinical-ready prognostic assay, treatment-selection model, or mechanistic validation package.
