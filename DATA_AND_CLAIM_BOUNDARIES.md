# Data and claim boundaries

This repository supports a discovery-stage, hypothesis-generating transcriptomic analysis. It does not provide clinical-deployment software or functional-validation evidence.

## Boundaries by analytical layer

| Layer | Supported interpretation | Not established |
| --- | --- | --- |
| Single-cell GlycoHigh state | Rank-balanced AUCell-defined malignant-hepatocyte state | Mechanistic driver status |
| MIF/SPP1 ligand analyses | Patient-aware expression association and inference-based candidate layers | Ligand secretion, receptor activation, direct immune suppression |
| Spatial transcriptomics | Spot-level co-enrichment after composition adjustment | Single-cell contact or direct signaling |
| Four-gene readout | Retrospective tissue-level association | Finalized clinical prognostic model |
| GSE235863 | Exploratory non-response association | Treatment-benefit prediction or decision support |
| inferCNV | Malignant-feature support | Perfect single-cell truth label |

## Locked Figure 2C rule

The active workflow uses `locked_source_data/single_cell/Fig2C_15391_tumor_hepatocyte_GlycoHigh_GlycoLow_FIXED.csv`. Two cells share the median AUCell score. Downstream scripts therefore read the locked rank-balanced assignment instead of recreating groups with a simple threshold comparison.
