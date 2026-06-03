# Revision reproducibility note

This folder documents revision-stage analyses that supplement the original public submission workflow. The goal is to make the manuscript-to-repository chain auditable without redistributing large public raw matrices or local intermediate RDS objects.

## Completed revision analyses

- Patient-internal rank-balanced grouping stability: `scripts/revision/patient_internal_median_sensitivity.py`; Supplementary Table S10.
- Discovery-cohort leave-four-out expression sensitivity: `scripts/revision/gse149614_leave_four_out_expression_python.py`; Supplementary Figure S32 and Supplementary Table S13.
- Independent GSE189903 HCC single-cell locked replication: `scripts/revision/gse189903_locked_replication_python.py`; Supplementary Figure S31 and Supplementary Table S12.
- Public GSE238264 Visium spatial-block empirical sensitivity: `scripts/revision/gse238264_spatial_block_permutation_python.py`; Supplementary Figure S33 and Supplementary Table S14.
- Stable GSE14520 Cox reconstruction with binary AFP coding: `scripts/revision/gse14520_stable_cox_python_validation.py`; Supplementary Table S11.

## Prepared but not claimed as completed

The following scripts are included only as templates or prepared analyses that require additional unavailable local inputs. They must not be described as completed unless outputs are generated, checked, and mapped in `MANUSCRIPT_TO_REPOSITORY_MAP.tsv`.

- `scripts/revision/spatial_block_permutation_inference.R`: RCTD-adjusted spatial-block empirical refit requiring the original RCTD spot-proportion table.
- `scripts/revision/gse14520_firth_cox_refit.R`: Firth Cox sensitivity requiring the GSE14520 source CSV from the frozen workflow.
- `scripts/revision/independent_single_cell_locked_replication.R`: generic Seurat-based independent single-cell replication template requiring a standardized Seurat object and malignant-cell flag.

## Output placement

Small completed outputs that are suitable for redistribution should be placed under:

```text
expected_outputs/revision/
```

Recommended files include:

```text
Supplementary_Table_S10.xlsx
Supplementary_Table_S11.xlsx
Supplementary_Table_S12.xlsx
Supplementary_Table_S13.xlsx
Supplementary_Table_S14.xlsx
GSE149614_leave_four_out_score_qc.csv
GSE189903_locked_replication_summary.csv
GSE238264_spatial_block_permutation_summary.csv
GSE14520_stable_cox_terms.csv
ANALYSIS_STATUS.tsv
PUBLIC_FROZEN_RELEASE_AUDIT.txt
```

## Boundary statement

Large public raw matrices, Visium images, Seurat objects, RCTD objects, participant-level clinical source tables, and local cache files are not redistributed. Users should download public data from GEO, TCGA/GDC, UCSC Xena, or the original public sources and provide local paths as documented in each script.

The four-gene tissue readout remains a retrospective discovery-stage tissue-level signal rather than a clinical assay. CellChat, NicheNet, LIANA, expression-level ligand-receptor compatibility, and spatial co-enrichment analyses remain transcriptomic inference or association support, not functional validation.
