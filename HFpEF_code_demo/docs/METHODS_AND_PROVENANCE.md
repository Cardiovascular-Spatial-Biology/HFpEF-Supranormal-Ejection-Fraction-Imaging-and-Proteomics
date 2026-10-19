# Methods and provenance

The workflow is based on the available analysis scripts and manuscript methods. `generate_demo.R` creates the demonstration data independently, without reading study files. Synthetic results have no biological interpretation.

## Provenance

| Component | Available source | Status |
|---|---|---|
| Paired exercise model | `1. dif_in_dif_limma_cystatin_PRODUCTION_1.R` | Analytical formula, pairing, contrast and moderation retained; explicit ID matching replaces positional assumptions. |
| Ternary response | `3_HFpEF_proteomics_ternary_plot.R` | Relative absolute-effect coordinates retained; Cartesian rendering replaces ggtern; labels/layout regenerated. |
| Tissue plots | `3. protein_tissue_plot_PRODUCTION_1.R` and `4_HFpEF_proteomics_tissue_enrichment.R` | Available annotation/significance logic consolidated; assay-ID bug corrected; plotting layout regenerated. Exact final visual provenance still unconfirmed. |
| Resting comparisons | Final manuscript's description of unadjusted comparisons | Reconstructed models, not a verified recovery of final S4 source or counts. |
| PCA / PLS-DA | Selected sections of `1_PCA_Protein.R` | Final S5 sources identified; unneeded exploratory branches removed; metadata and X-variance calculation corrected. |
| CAMERA | Available `B4h_pathway_map_standalone.R` plus final written methods | Reconstructed implementation. Final author script, reference releases and figure output have not been matched. |
| Pathway enrichment | Available `2_HFpEF_proteomics_reactome_pathways.R`, older Enrichr calls, final written methods | Reconstructed Fig 4d. Offline fictional-set demonstration plus custom-background Enrichr export/live modes; no claimed recovery of final study results. |
| Medication sensitivity | Final manuscript, response text and Table S14 caption | Reconstructed within-group paired interaction models. Prevalence-centred exercise contrast and retention definitions are explicit choices pending final-author confirmation. |

## Paired exercise response and overlap

For each group separately, use log2 RFU, design `~0 + Condition + Age + Sex + BMI + BaselineCystatin`, a participant block, `limma::duplicateCorrelation`, `lmFit`, the Stress minus Rest contrast, and `eBayes` with default moderation. Baseline cystatin is the log2 resting signal of the designated assay. No covariate is silently omitted if unavailable. BH correction is applied over assays within each group.

Exercise-responsive sets require FDR <0.05 and `abs(log2FC) > log2(1.25)`. This exact cutoff (approximately 0.321928) retains the original set-selection rule, rather than rounding it to 0.32. The Venn diagram combines increases and decreases. An assay significant in one group and not another is **not by itself evidence of a significant between-group interaction**.

Ternary coordinates are `(|log2FC_g| + 0.01) / sum_g(|log2FC_g| + 0.01)`. They represent relative magnitudes and do not encode direction. The axes have Control on the left, higherEF at the top and lowerEF on the right. Point sizes show the maximum absolute log2FC. Gene labels are selected by extreme position for a readable demonstration, not manually selected study targets.

## Tissue annotation and plotting

Pipe-delimited UniProt mappings are expanded while retaining the original `SeqId`. The prior assignment of a data-frame row number to `SeqId` was removed. Repeated assay–tissue pairs are plotted once; the exported mapping records this choice. Multiple tissues for one assay can still produce multiple points.

The supplied nTPM summary is parsed. A single listed tissue is accepted; otherwise the largest listed value must be at least five times the second largest, with a positive maximum. This retains the available code's summary-field rule. **It does not establish that the selected tissue is five times every tissue in a complete atlas.** A complete atlas-based definition would require separately verified reference data and a documented change to the annotation method.

For each of the two HFpEF panels, coloured assays require FDR <0.05 in that panel, without a fold-change cutoff. Shared categories additionally require the same direction of change across the significant groups. Thus tissue colours intentionally differ from the overlap/ternary threshold. Y is `-log10(FDR)`, size is absolute log2FC, and shape is increase/decrease. Coordinates and label placement are regenerated rather than pixel-identical to the manuscript.

## Resting comparisons

Using resting log2 RFU, fit `~Group` without covariates, separately for pooled HFpEF minus Control, and higherEF minus lowerEF among HFpEF participants. Apply empirical Bayes moderation and BH adjustment independently for each comparison; significance is FDR <0.05. These choices follow the written methods/captions. No claim is made that the original final S4 plotting conventions or numerical outputs have been recovered.

## PCA and PLS-DA

Retain `log2(RFU + 1)` from the original S5 code, which differs from the plain log2 transform in the limma/CAMERA modules. Resting PCA removes zero-variance assays and centres/scales variables. PLS-DA uses all assays and all samples, two components, unit-variance scaling, and either condition or phenotype labels. No differential-protein or top-variable selection precedes these models.

Sample labels are matched directly from `sampledata.csv`; no `Sample.x` or `Condition.x` from an old R session is required. Both PLS-DA panels label X-variance using `fit$prop_expl_var$X`. The previous S5c custom score-variance ratio did not represent the model's explained X-variance. Correcting it changes the axis percentages, not the fitted scores.

These are exploratory in-sample displays. PLS-DA knows the labels during fitting; separation alone is not independent evidence of predictive performance. There is no covariate adjustment or within-subject multilevel decomposition in these visualisations. The ellipses display score distributions, not tests of group differences. If prediction is evaluated later, split by participant so paired samples remain in the same fold. See the [mixOmics performance guidance](https://mixomics.org/performance-assessment-and-parameter-tuning/) and [multilevel guidance](https://mixomics.org/methods/multilevel/).

## CAMERA reconstruction

1. Restrict participants to the two HFpEF groups.
2. Keep unambiguous gene labels. Select the assay with highest mean **raw resting RFU among HFpEF participants**, breaking exact ties by assay ID. This population/scale choice is explicit and still needs confirmation against the final author's implementation.
3. Calculate each participant's `log2(stress) - log2(rest)` directly. Do not difference fitted values.
4. For each of the seven specified clinical/imaging variables, use complete cases and the design `~x + Age + Sex + BMI`. There is no cystatin covariate in this pathway design, matching the written methods.
5. Use supplied, versioned pathway membership at levels 2–3; retain 10–500 measured genes. The code consumes the provided hierarchy level rather than reconstructing a Reactome hierarchy itself. Constant responses within a variable's complete cases are removed, and measured set sizes are recorded.
6. Call `limma::camera` with **`inter.gene.cor=0.01` explicitly**, retaining the agreed primary setting. BH adjustment is within variable. An additional global BH column is descriptive and does not determine the figure.
7. If a variable has FDR-significant sets, display those; otherwise display at most ten nominal p <0.05 sets. Family membership is an explicit input, not inferred from substrings. Ribbon width counts contributing pathways. Because sets can overlap, width is not independent evidence. Direction is positive if at least 75% of contributing sets are Up, negative if at most 25% are Up, otherwise mixed. Pale ribbons are nominal-only; dark ribbons contain at least one FDR-significant set.

The new module includes LA end-systolic volume index change instead of the obsolete resting NT-proBNP predictor and performs one-assay-per-gene selection. It does **not** reproduce or enforce the manuscript's number of pathways or significant associations. Exact Reactome release, hierarchy export, final phenotype version, family assignments and final source are still required to validate study reproduction. The supplied example has fictional pathways and labels, so its associations have no biological interpretation.

## Fig 4d and Table S14

Fig 4d uses all unique unambiguous measured genes as the background, with significant assays mapped to unique query genes. The default offline demonstration uses a hypergeometric upper-tail test and BH over all tested fictional sets, including zero overlaps. It is not presented as Enrichr. For study use, cached or live Enrichr modes preserve the service's adjusted p-values and require the specified `Reactome_Pathways_2024` library and measured background. Plot exclusions and priority grouping follow the available plotting source, with the ambiguous “Control” category relabelled as “Control-containing (other)”. Full results are retained before display filtering.

Table S14 fits each HFpEF group and drug separately with the paired model, age, sex, BMI, baseline cystatin, medication main effect and medication-by-exercise interaction. Medication is centred at its observed prevalence among complete cases; the exercise coefficient therefore estimates a prevalence-weighted average response. This avoids inadvertently comparing the primary population average with the non-user-only response. The precise final author's contrast was not available, so this centering choice is documented rather than claimed to match it. The interaction tests the difference in exercise response between users and non-users. BH correction uses all assays within each group/drug/contrast.

The primary set uses both FDR and fold change. The displayed retention count uses adjusted FDR alone, consistent with the response text's stated retention criterion, while complete CSVs also report retention using both thresholds, direction and new hits. Missing medication values lead to matched complete-case baseline and adjusted fits, with denominators reported. Unestimable strata are labelled explicitly. Full details and optional Enrichr instructions are in [FIG4D_AND_TABLES14.md](FIG4D_AND_TABLES14.md).

## Scope

The old CCA, unused quadrant plots, tissue-count comparisons, UMAP, t-SNE, clustering, selected-protein exploratory PCAs and unused six-class PLS-DA were not carried forward. Clinical-imaging statistics, pressure-volume modelling and digital-twin simulations remain outside this package.

## Input validation and reproducibility

Inputs use explicit participant, sample and assay keys. Validation rejects duplicate IDs, inconsistent matrices, invalid RFU, missing adjustment variables and non-identifiable designs. Each run records software versions and input checksums without filesystem locations or a full session dump.

The distributed data and reference annotations are generated de novo. Study matrices, clinical records and cached study models are not included.
