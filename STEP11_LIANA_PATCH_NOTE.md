# STEP11_LIANA_PATCH_NOTE

## Purpose

This Step 11 patch adds a minimal LIANA directional-concordance analysis requested during final cross-check. It resolves a Methods/Results/Figure legend/Supplementary consistency gap without changing the manuscript positioning.

## Analysis scope

The analysis compares predefined MIF- and SPP1-related ligand-receptor axes in the following restricted sender-receiver framework:

- Sender groups: GlycoHigh tumor-derived hepatocytes and GlycoLow tumor-derived hepatocytes
- Receiver groups: tumor T/NK-lineage cells and tumor myeloid cells
- Candidate axes: MIF-CD74, MIF-CXCR4, MIF-CD44, SPP1-CD44, SPP1-ITGAV, SPP1-ITGB1, SPP1-ITGA4, SPP1-ITGB5, and SPP1-ITGB6

## Key output

Among 18 predefined axis-receiver combinations, 14 showed stronger aggregate ranks for GlycoHigh hepatocytes than for GlycoLow hepatocytes, four showed equal ranks, and none favored GlycoLow hepatocytes.

## Interpretation boundary

LIANA was used only as a cross-method directional-concordance check. It must not be interpreted as:

- functional ligand-receptor validation
- mechanistic proof
- ligand secretion evidence
- receptor activation evidence
- spatial proximity evidence
- direct signaling evidence
- treatment-response evidence
- clinical decision support

## Manuscript wording

Use: candidate / inference-based / directional concordance / cross-method support.

Avoid: validated communication / functional validation / proven signaling / drive / mediate / mechanistic validation.

## Related files

- `15_LIANA_directional_concordance.R`
- `results/LIANA_directional_concordance/LIANA_group_counts.csv`
- `results/LIANA_directional_concordance/Supplementary_Table_S4_LIANA_full_results.tsv`
- `results/LIANA_directional_concordance/Supplementary_Table_S4_LIANA_aggregate_results.tsv`
- `results/LIANA_directional_concordance/Supplementary_Table_S4_LIANA_candidate_axes.tsv`
- `results/LIANA_directional_concordance/LIANA_direction_counts.tsv`
- `results/LIANA_directional_concordance/Figure4D_LIANA_directional_concordance.pdf`
- `results/LIANA_directional_concordance/Figure4D_LIANA_directional_concordance.png`
- `results/LIANA_directional_concordance/sessionInfo_LIANA.txt`
