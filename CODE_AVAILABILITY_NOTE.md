# Code availability note

Analysis scripts, compact locked source tables, revision-stage expected outputs, and reproducibility documentation are available in this repository.

## Repository

Repository URL:

`https://github.com/17698988000/HCC-glycolysis-immune-scRNAseq`

Submission revision final release:

`v1.0.2-submission-revision-final`

Release commit:

The immutable release commit is recorded on the GitHub release/tag page for `v1.0.2-submission-revision-final`.

## Included

- Active analysis scripts under `scripts/`.
- Revision-stage analysis and figure-redraw scripts under `scripts/revision/`.
- Compact revision-stage expected outputs under `expected_outputs/revision/`.
- Authoritative submission-locked compact source tables under `locked_source_data/`.
- Current author environment snapshot and reproducibility files under `reproducibility/`.
- Revision reproducibility explanation under `docs/REVISION_REPRODUCIBILITY_NOTE.md`.
- Recommended execution order in `RUN_ORDER.md`.
- Manuscript-to-repository mapping in `MANUSCRIPT_TO_REPOSITORY_MAP.tsv`.
- Historical or unused materials retained separately under `archive_not_used/`.

## Excluded

Large downloaded public datasets, Visium image files, TCGA/Xena downloads, intermediate RDS objects, local cache files, and generated manuscript figures are not committed.

Public dataset accessions are documented in the manuscript and README. Manuscript figures are submitted separately with the journal package.

## Interpretation boundary

This repository supports reproducibility and auditability of the submitted manuscript.

It does not provide a clinical assay, validated prognostic classifier, treatment-response predictor, or mechanistically validated ligand-receptor model.

## Environment snapshot note

The files in `reproducibility/` capture the author's current environment at repository-cleanup time. They are provided transparently and may differ from the original analysis environment reported in the manuscript. The manuscript remains the authoritative record of the analysis environment used for the study.
