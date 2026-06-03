# Recommended run order

Run scripts from the repository root. The repository contains compact locked source tables and analysis scripts; large public-data downloads and local intermediate RDS files are intentionally excluded.

## Main submission workflow

| Step | Script | Role |
| ---: | --- | --- |
| 1 | `scripts/01_QC_clustering.R` | Upstream single-cell QC, integration, clustering, and annotation |
| 2 | `scripts/02_glycolysis_scoring.R` | Validate and export the authoritative rank-balanced Figure 2C GlycoHigh/GlycoLow assignment |
| 3 | `scripts/06_inferCNV_malignant.R` | Malignant-feature marker and inferCNV support |
| 4 | `scripts/16_patient_level_ligand_effects.R` | Patient-aware MIF/SPP1/PTGES summaries aligned to Figure 3 and Supplementary Table S3 |
| 5 | `scripts/03_cellchat_analysis.R` | CellChat inference-based communication support |
| 6 | `scripts/15_NicheNet_analysis.R` | NicheNet transcriptional-response sensitivity analyses |
| 7 | `scripts/18_LIANA_directional_concordance.R` | LIANA directional-concordance check for predefined MIF/SPP1 axes |
| 8 | `scripts/11_spatial_transcriptomics.R` | Spatial transcriptomic and RCTD-related support |
| 9 | `scripts/07_LASSO_risk_score.R` | TCGA four-gene tissue-level readout and repeated-CV stability analysis |
| 10 | `scripts/10_GSE14520_validation.R` | External fixed-coefficient survival evaluation |
| 11 | `scripts/12_GSE235863_exploratory_association.R` | Exploratory anti-PD-1 plus lenvatinib association only |
| 12 | `scripts/08_glycolysis_gradient.R` | Locked Supplementary Figure S15 redraw from compact source table |
| 13 | `scripts/09_TF_activity.R` | Auxiliary TF-activity analysis |
| 14 | `scripts/14_OXPHOS_metabolic_specificity.R` and `scripts/17_partial_correlation_metabolic_specificity.R` | Auxiliary metabolic-specificity analyses |

## Revision sensitivity workflow

| Step | Script | Role |
| ---: | --- | --- |
| 15 | `scripts/revision/patient_internal_median_sensitivity.py` | Supplementary Table S10; patient-internal rank-balanced grouping stability |
| 16 | `scripts/revision/gse14520_stable_cox_python_validation.py` | Supplementary Table S11; stable GSE14520 Cox reconstruction with binary AFP coding |
| 17 | `scripts/revision/gse189903_locked_replication_python.py` | Supplementary Figure S31 and Supplementary Table S12; independent GSE189903 locked replication |
| 18 | `scripts/revision/gse149614_leave_four_out_expression_python.py` | Supplementary Figure S32 and Supplementary Table S13; discovery-cohort leave-four-out expression sensitivity |
| 19 | `scripts/revision/gse238264_spatial_block_permutation_python.py` | Supplementary Figure S33 and Supplementary Table S14; public-Visium spatial-block empirical sensitivity |

## Prepared templates and diagnostic utilities

The following scripts are included for transparency or future checking, but should not be treated as completed manuscript-supporting analyses unless their outputs are explicitly mapped:

- `scripts/revision/independent_single_cell_locked_replication.R`
- `scripts/revision/leave_four_out_patient_internal_expression_sensitivity.R`
- `scripts/revision/spatial_block_permutation_inference.R`
- `scripts/revision/gse14520_firth_cox_refit.R`
- `scripts/revision/download_geo_revision_inputs.py`
- `scripts/revision/generate_gse149614_leave_four_out_figure.py`
- `scripts/revision/generate_gse189903_replication_figure.py`
- `scripts/revision/generate_gse238264_spatial_block_figure.py`
- `scripts/revision/generate_revised_figure6.py`

## Notes

- `scripts/utils/locked_fig2c_groups.R` is the authoritative downstream grouping helper.
- Two cells share the Figure 2C median AUCell score. Do not recreate groups using a simple `> median` comparison.
- Manuscript figures are submitted separately and are not required in the code repository.
- Revision-stage compact expected outputs are stored in `expected_outputs/revision/`.
