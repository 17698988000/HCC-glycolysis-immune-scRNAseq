# Submission release checklist

## Repository structure

- [x] Active scripts moved under `scripts/`
- [x] Revision-era scripts omitted from the public release and documented in `archive_not_used/ARCHIVED_ITEMS_MANIFEST.csv`
- [x] Compact authoritative source tables added under `locked_source_data/`
- [x] MIT License and citation metadata added
- [x] Machine-specific `D:/scRNA_project` paths removed from active scripts
- [x] Figure 2C downstream grouping aligned to locked rank-balanced assignment
- [x] Supplementary Figure S15 redraw aligned to 15,391-cell equal-cell-count-bin source table
- [x] Figure 3 ligand QC targets aligned to verified final values

## Before publishing GitHub release

- [ ] Review static-validation report supplied with the candidate package
- [ ] Sync candidate repository to the `submission-cleanup` branch
- [ ] Confirm GitHub homepage rendering
- [ ] Create pull request or merge `submission-cleanup` into `main`
- [ ] Create release tag: `v1.0.0-submission`
- [ ] Copy final commit hash into manuscript Code availability statement
- [ ] Confirm Supplementary Figure S15 vector PDF remains in the journal submission package

## Environment note

- [ ] Confirm whether an original R 4.4.2 session snapshot is available. The committed snapshot currently reflects the author's local cleanup-time R 4.5.2 environment and is labeled accordingly.
