# Code availability note

This file documents the code-availability statement for the current Step-11 cross-checked manuscript package.

## Repository URL

All R analysis scripts used in this study are available at:

```text
https://github.com/17698988000/HCC-glycolysis-immune-scRNAseq
```

## Recommended manuscript wording

All R analysis scripts used in this study are available at: https://github.com/17698988000/HCC-glycolysis-immune-scRNAseq.

The revised script set follows the Figure 1-7 manuscript structure and includes single-cell preprocessing, AUCell glycolysis scoring, malignant-feature and inferCNV/CNV support, patient-level ligand summaries, CellChat analysis, NicheNet sensitivity analysis, focused ligand-receptor expression compatibility checks, LIANA directional-concordance analysis for predefined MIF/SPP1 axes, RCTD-based spatial deconvolution and composition-adjusted spatial regression, TCGA-LIHC and GSE14520 survival analyses, proportional-hazards diagnostics, time-dependent ROC analysis, calibration analysis, bootstrap optimism correction, benchmarking, immune-context analyses, and the GSE235863 exploratory non-response association analysis.

Analyses were performed in R version 4.4.2; key analytical packages and statistical frameworks are described in the Methods.

## Interpretation boundary

This repository supports a discovery-stage and hypothesis-generating transcriptomic study. It should not be interpreted as a clinical-ready prognostic assay, treatment-selection model, or mechanistic validation package.

The TPI1/ENO1/LDHA/SLC2A1 score is a tissue-level readout of GlycoHigh biology, not an implemented clinical assay. SPP1 and MIF analyses provide candidate, inference-based, complementary communication-layer support only. GSE235863 is used only for exploratory association analysis and not for treatment-benefit inference, treatment selection, or clinical decision-making.
