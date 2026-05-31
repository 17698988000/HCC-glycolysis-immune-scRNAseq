# Recommended run order

Run scripts from the repository root. The repository contains compact locked source tables and analysis scripts; large public-data downloads and local intermediate RDS files are intentionally excluded.

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

## Notes

- `scripts/utils/locked_fig2c_groups.R` is the authoritative downstream grouping helper.
- Two cells share the Figure 2C median AUCell score. Do not recreate groups using a simple `> median` comparison.
- Manuscript figures are submitted separately and are not required in the code repository.
