# Submission release checklist

## Repository structure

- [x] Active scripts organized under `scripts/`.
- [x] Revision-stage scripts included under `scripts/revision/`.
- [x] Compact revision-stage expected outputs included under `expected_outputs/revision/`.
- [x] Compact authoritative source tables retained under `locked_source_data/`.
- [x] README homepage renders correctly.
- [x] `RUN_ORDER.md` includes main and revision-stage execution order.
- [x] `MANUSCRIPT_TO_REPOSITORY_MAP.tsv` maps manuscript figures/tables to scripts and compact outputs.
- [x] `docs/REVISION_REPRODUCIBILITY_NOTE.md` explains completed revision analyses and prepared-but-not-claimed templates.
- [x] `.gitignore` excludes large raw data, RDS objects, archives, cache files, and local environments.
- [x] MIT license and citation metadata included.
- [x] Machine-specific local paths are avoided in public-facing documentation.

## Submission-locked analysis boundaries

- [x] Figure 2C downstream grouping aligned to locked rank-balanced assignment.
- [x] Four-gene tissue-level readout documented as discovery-stage tissue translation, not a clinical assay.
- [x] MIF/SPP1, CellChat, NicheNet, LIANA, and spatial analyses described as inference or compatibility evidence, not functional validation.
- [x] Large public raw matrices, Visium image files, TCGA/Xena downloads, local cache files, and large intermediate RDS objects are not redistributed.

## Before publishing final GitHub release

- [ ] Confirm GitHub homepage rendering after the latest README update.
- [ ] Confirm `RUN_ORDER.md`, `MANUSCRIPT_TO_REPOSITORY_MAP.tsv`, and `docs/REVISION_REPRODUCIBILITY_NOTE.md` render correctly.
- [ ] Confirm `scripts/revision/` and `expected_outputs/revision/` are visible in the public repository.
- [ ] Confirm no `__pycache__`, `.pyc`, `.rds`, `.RDS`, `.zip`, `.tar`, `.gz`, `.h5`, `.h5ad`, `.mtx`, or local raw-data files are present in the public release.
- [ ] Create or update release tag `v1.0.0-submission` after all final commits.
- [ ] Copy the final release link and final commit hash into the manuscript Code availability statement.
- [ ] Confirm submitted manuscript figures remain in the journal submission package and are not required to be rebuilt from the GitHub repository.

## Environment note

- [ ] Confirm the final committed session/environment notes match the manuscript and repository claim boundaries.
- [ ] Confirm `requirements-revision.txt` is present for Python revision-stage scripts.
