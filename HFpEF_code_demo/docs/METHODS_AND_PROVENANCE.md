# Methods

The paired limma model retains the available study implementation. 

## Exercise response and overlap

Within each group, fit log2 RFU with `~0 + Condition + Age + Sex + BMI + BaselineCystatin`, participant blocking, `limma::duplicateCorrelation` and empirical Bayes moderation. Test Stress minus Rest; apply BH correction over all assays within each group. Baseline cystatin is the log2 resting signal of the designated assay.

Venn and ternary categories use FDR <0.05 and `abs(log2FC) > log2(1.25)`. Ternary coordinates are `(|log2FC_g| + 0.01) / sum_g(|log2FC_g| + 0.01)` and show magnitude, not direction. 

## Tissue plots

Expand UniProt mappings while retaining assay IDs. Accept a single listed tissue or a maximum nTPM at least five times the second-largest listed value. This rule applies to the supplied summary, not a complete atlas. Colours use FDR <0.05 without a fold-change cutoff; shared categories require matching response direction.

## Resting comparisons and PCA/PLS-DA

Resting comparisons use log2 RFU, unadjusted `~Group` models and BH correction per comparison: pooled HFpEF versus Control, and higherEF versus lowerEF.

PCA and PLS-DA use centred, scaled `log2(RFU + 1)`. PCA uses resting samples and removes zero-variance assays. PLS-DA uses all samples and assays with two components; axis labels report X-variance. 

## CAMERA

Use HFpEF participants and one assay per unambiguous gene, selected by highest mean raw resting RFU. Calculate `log2(stress) - log2(rest)` for each participant. For each variable in `config.R`, use complete cases and `~x + Age + Sex + BMI`.

Test supplied level 2-3 pathways containing 10-500 measured genes with `inter.gene.cor = 0.01`. BH correction is within variable. Display FDR <0.05 pathways, or up to ten nominal p <0.05 pathways if none pass FDR. Dark ribbons contain an FDR-significant pathway; pale ribbons are nominal-only. Width counts pathways, which may overlap. Direction is positive when at least 75% are Up, negative when at most 25% are Up, otherwise mixed.

## Medication sensitivity

See [Table S14](TABLES14.md) for the paired interaction model and retention definition.
