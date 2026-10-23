# Table S14: medication sensitivity

`R/08_medication_sensitivity.R` tests ACEi, ARB, MRA, BetaBlocker, SGLT2i and LoopDiuretic separately within each HFpEF group. Medication fields use 1=yes, 0=no and NA=unknown.

Baseline and adjusted models use the same complete-case participants, log2 RFU, participant blocking, `duplicateCorrelation` and empirical Bayes moderation. The adjusted design is:

```r
DrugCentered <- Drug - mean(Drug)
~ Condition * DrugCentered + Age + Sex + BMI + BaselineCystatin
```

The exercise coefficient estimates the response at the observed medication prevalence. The interaction estimates the user-minus-non-user difference in exercise response. BH correction is over all assays within each group, drug and contrast. This reconstructed contrast has not been matched to the final author model.

The baseline set requires FDR <0.05 and |log2FC| >log2(1.25). Retention counts baseline assays with adjusted FDR <0.05 on the same participants. CSV outputs also report fold-change retention, direction, interactions and sample counts. Fewer than two users or non-users, or a rank-deficient model, is marked not estimable.
