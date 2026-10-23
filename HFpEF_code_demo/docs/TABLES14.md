# Table S14: methods and usage

This analysis runs as part of the default seven-step demonstration. Results are calculated from the supplied synthetic inputs. The final medication-analysis source was unavailable, so numerical agreement with the final manuscript has not been established.

## Table S14

`R/08_medication_sensitivity.R` analyses the six classes described in the supplement: ACE inhibitor, ARB, MRA, beta blocker, SGLT2 inhibitor and loop diuretic. The corresponding numeric fields are `ACEi`, `ARB`, `MRA`, `BetaBlocker`, `SGLT2i`, `LoopDiuretic`: 1=use, 0=no use, NA=unknown. Classes are entered one at a time, within each HFpEF group.

For each class, retain participants with known medication status, fit both baseline and adjusted models on those same participants, and estimate the repeated-measure correlation with `duplicateCorrelation`. Fit log2 RFU with participant blocking and empirical Bayes moderation. The adjusted design is:

```r
DrugCentered <- Drug - mean(Drug)
~ Condition * DrugCentered + Age + Sex + BMI + BaselineCystatin
```

`ConditionStress` is the exercise effect at the group’s observed medication prevalence. In an equivalent uncentred model, it is `beta_Stress + prevalence * beta_Stress:Drug`. `ConditionStress:DrugCentered` is the user-minus-non-user difference in exercise response. Both contrasts are tested over all assays, with BH correction within group, drug and contrast. An additional pooled interaction FDR is exported but does not drive the main table.

This centering is a documented reconstruction choice. If the final author used a non-user contrast, equal weighting, or another standardisation, numerical results will differ. Neither the manuscript's reported retained counts nor its lack of significant interactions is enforced in this code. A non-significant interaction is not proof of no medication effect or no confounding.

The baseline set requires FDR <0.05 and |log2FC| >log2(1.25). The main table reports how many baseline assays retain FDR <0.05 after medication adjustment, divided by the baseline set size on the same complete cases. Detailed CSVs additionally show retention with both FDR and fold change, retention with the same direction, and newly significant assays. With missing medication values, the denominator may differ from the full primary-model count; sample counts and missingness are explicit.

Strata with fewer than two users or non-users, or a rank-deficient design, are marked **not estimable** rather than assigned zero interactions. Complete outputs include the display CSV, a PDF/PNG table, all assay-level adjusted responses and interactions, sample/model diagnostics and a caption defining the estimand.

## Synthetic identifiers

`SYN` stands for **synthetic**. `SYN_P001` means fictional participant 001, not a disease code or a pseudonym linked to a real person. Keeping this prefix helps prevent accidental confusion between demonstration and study files. The six medication fields are generated from independent random numbers; no patient treatment records were used.
