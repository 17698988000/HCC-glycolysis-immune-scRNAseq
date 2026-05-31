# HCC GlycoHigh malignant-hepatocyte state

Analysis code and locked source tables for a discovery-stage hepatocellular carcinoma transcriptomic study integrating single-cell, spatial, and retrospective bulk-cohort analyses.

## Study scope

This repository supports a **hypothesis-generating** analysis package. It is not a clinical-ready prognostic assay, treatment-selection model, or mechanistic-validation package.

Key interpretation boundaries:

- The TPI1/ENO1/LDHA/SLC2A1 score is a tissue-level readout of GlycoHigh biology, not an implemented clinical assay.
- SPP1 and MIF are candidate, inferred, complementary communication or transcriptional-response layers.
- CellChat, NicheNet, expression-level ligand-receptor compatibility, and LIANA outputs are transcriptomic inference support, not functional validation.
- GSE235863 is used only for exploratory association analysis in a small, class-imbalanced cohort.
- inferCNV is used as malignant-feature support, not as a perfect single-cell truth label.

## Submission-locked conventions

The authoritative Figure 2C assignment is stored in:

```text
locked_source_data/single_cell/Fig2C_15391_tumor_hepatocyte_GlycoHigh_GlycoLow_FIXED.csv
```

Locked single-cell scope:

| Item | Locked value |
| --- | ---: |
| Tumor-derived hepatocytes | 15,391 |
| Median AUCell glycolysis score | 0.212617779598525 (displayed as 0.213) |
| Cells tied at the median | 2 |
| GlycoHigh | 7,695 |
| GlycoLow | 7,696 |
| Grouping rule | Rank-balanced split using the locked Figure 2C assignment table |

The locked TCGA-derived tissue-level readout is:

```text
score = 0.3041908 * TPI1 +
        0.9639654 * ENO1 +
        1.3404374 * LDHA +
        0.2424239 * SLC2A1
```

The lambda.min model retained four genes. The conservative lambda.1se model was null and is treated as an important limitation.

## Repository structure

```text
scripts/                 active analysis scripts
scripts/utils/           shared locked-assignment helper
locked_source_data/      compact authoritative source tables used for the submission release
reproducibility/         current author environment snapshot and refresh script
docs/                    method-specific notes
archive_not_used/        manifest of revision-era materials omitted from the public release
```

Generated figures and large intermediate objects are intentionally not committed. Manuscript figures are submitted separately; scripts and compact locked source tables are retained here for reproducibility and provenance.

## Public data resources

| Dataset | Manuscript role | Boundary |
| --- | --- | --- |
| GSE149614 | Single-cell discovery | Eight paired HCC patients after filtering |
| GSE238264 | Spatial support | Four Visium sections; spot-level co-enrichment only |
| TCGA-LIHC | Bulk discovery and clinical association | Retrospective tissue-level association |
| GSE14520 | External survival evaluation | Fixed TCGA-derived coefficients |
| GSE235863 | Exploratory treatment-associated analysis | Small, class-imbalanced; hypothesis-generating only |

## Reproducing the workflow

Run scripts from the repository root. Large public-data downloads and intermediate RDS objects are not committed. Set local paths with environment variables where needed, for example:

```r
Sys.setenv(SEURAT_FINAL_RDS = "path/to/seurat_final.rds")
Sys.setenv(PROJECT_DIR = "path/to/repository")
```

See [`RUN_ORDER.md`](RUN_ORDER.md) for the recommended execution sequence and [`MANUSCRIPT_TO_REPOSITORY_MAP.tsv`](MANUSCRIPT_TO_REPOSITORY_MAP.tsv) for the figure-to-script map.

## License

Code is released under the [MIT License](LICENSE).
