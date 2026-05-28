# Recommended run order for the current manuscript

This order follows the Step-8-final manuscript and the Step-9 supplementary-repository requirement. It is not intended to create new claims or new manuscript figures without a documented QC pass.

1. `01_QC_clustering.R`  
   Build the upstream single-cell object if reprocessing from raw public data.

2. `02_glycolysis_scoring.R`  
   Confirm locked tumor-derived hepatocyte counts and GlycoHigh/GlycoLow grouping.

3. `06_inferCNV_malignant.R`  
   Generate marker-score source/QC outputs and optional inferCNV/CNV support when local inferCNV dependencies are present.

4. `16_patient_level_ligand_effects.R`  
   Generate patient-level ligand robustness summaries for SPP1, MIF, and PTGES.

5. `03_cellchat_analysis.R` and `15_NicheNet_analysis.R`  
   Generate or guard inference-based communication and ligand-activity outputs. Do not interpret as functional validation.

6. `14_OXPHOS_metabolic_specificity.R` and `17_partial_correlation_metabolic_specificity.R`  
   Generate metabolic-specificity source/QC outputs.

7. `11_spatial_transcriptomics.R` plus retained scripts in `revision_figure_restore/`  
   Generate Visium/RCTD spatial support and retained spatial supplementary/model-diagnostic visualizations.

8. `07_LASSO_risk_score.R` and `10_GSE14520_validation.R`  
   Generate TCGA-derived four-gene readout and external score-level survival evaluation. Continuous Cox modelling remains the primary survival framework.

9. `12_GSE235863_exploratory_association.R`  
   Generate descriptive GSE235863 source/QC outputs only. This script does not run ROC/PR/DCA/ML.

10. `13_revision_analyses.R`  
   Optional status/index check only. Do not use it to create new claims.

Archived files under `archive_not_used/` are not part of the current manuscript run order.
