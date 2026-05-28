# Data and claim boundaries

## Dataset roles

| Dataset | Current role | Boundary |
|---|---|---|
| GSE149614 | single-cell discovery layer | Defines GlycoHigh malignant-hepatocyte state; retrospective transcriptomic evidence only. |
| GSE238264 | spatial transcriptomics support | Spot-level tissue co-enrichment after composition adjustment; not same-cell expression or direct ligand-receptor contact. |
| TCGA-LIHC | bulk discovery and clinical association | Four-gene tissue-level readout construction and score-level survival association; not a clinical-ready assay. |
| GSE14520 | external score-level survival evaluation | External evaluation using TCGA-derived coefficients; not clinical validation. |
| GSE235863 | exploratory anti-PD-1 plus lenvatinib association | Underpowered (n = 15; non-responders = 4); descriptive/hypothesis-generating only. |
| GSE125449 | archived / not used in current manuscript | Removed from current dataset census; retained only in `archive_not_used/` for historical traceability. |

## Claim boundaries

- Use: discovery-stage, hypothesis-generating, tissue-level readout, candidate/inferred/complementary communication layers, CNV/malignancy support.
- Avoid: validated predictor, clinical assay, clinical classifier, treatment-selection model, functional validation, causal mechanism, same-cell co-expression, clinical implementation.

## GSE235863 note

The current manuscript removed ROC/PR/performance-oriented GSE235863 analysis. The active script is limited to score distribution, median-stratified descriptive summary, Wilcoxon comparison, and Fisher exact summary.
