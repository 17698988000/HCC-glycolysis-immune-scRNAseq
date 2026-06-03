# Revision reproducibility note

This folder documents revision-stage reproducibility materials added for the CSBJ revision.

## Completed revision analyses

The following revision-stage analyses are treated as completed and mapped in the repository:

- Patient-internal GlycoHigh/GlycoLow grouping stability: Supplementary Table S10.
- Stable GSE14520 Cox reconstruction with binary AFP coding: Supplementary Table S11.
- Independent GSE189903 locked single-cell replication: Supplementary Figure S31 and Supplementary Table S12.
- GSE149614 leave-four-out expression sensitivity excluding TPI1, ENO1, LDHA, and SLC2A1: Supplementary Figure S32 and Supplementary Table S13.
- GSE238264 public-Visium spatial-block empirical sensitivity: Supplementary Figure S33 and Supplementary Table S14.

## Prepared but not claimed as completed

The following scripts may be included as templates or diagnostic utilities, but they should not be interpreted as completed manuscript-supporting analyses unless their outputs are explicitly mapped:

- Generic Seurat-based independent single-cell replication template.
- RCTD-adjusted spatial-block empirical refit template.
- GSE14520 Firth Cox sensitivity template.

## Data redistribution boundary

Large public raw matrices, Visium image files, downloaded GEO archives, TCGA/Xena files, large RDS objects, and local cache files are not redistributed in this repository.

Users should download public data from GEO, TCGA/GDC, UCSC Xena, or other original sources and provide local paths as documented in the scripts.

## Interpretation boundary

These revision analyses support reproducibility and auditability of the submitted manuscript. They do not convert the four-gene readout into a clinical assay, validated prognostic classifier, treatment-response predictor, or mechanistically validated ligand-receptor model.
