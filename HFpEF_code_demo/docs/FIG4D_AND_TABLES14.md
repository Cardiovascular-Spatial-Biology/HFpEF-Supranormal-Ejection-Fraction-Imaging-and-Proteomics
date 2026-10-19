# Fig 4d and Table S14: methods and usage

Both analyses run as part of the default eight-step demonstration. Results are calculated from the supplied synthetic inputs. The final study enrichment exports and medication-analysis source were unavailable, so numerical agreement with the final manuscript has not been established.

## Fig 4d

`R/07_pathway_enrichment.R` builds a unique measured-gene background from feature annotations. Each group’s query is the union of gene symbols represented by assays with FDR <0.05 and absolute log2FC >log2(1.25), combining increases and decreases. Multi-gene annotations are excluded explicitly; repeated assays do not count repeatedly as genes. All query/background lists and their hashes are exported.

Three modes are available through `config.R`:

| Mode | What is run | Network |
|---|---|---|
| `demo_local` (default) | Hypergeometric upper-tail enrichment of fictional `DEMO_` sets against fictional measured genes; BH over all nonempty background-mapped sets, including zero-overlap sets | None |
| `enrichr_exports` | Load original Enrichr group exports, preserving the service's raw and adjusted p-values | None |
| `enrichr_live` | Call the named Enrichr library with an explicit measured-gene background and cache results | Only with `--allow-enrichr-network` |

The default is **not an Enrichr reproduction**. It keeps the synthetic demonstration runnable without sending fictional gene labels to an external service. It does not impose CAMERA’s 10–500 gene size or hierarchy filters on Fig 4d. The service modes specify `Reactome_Pathways_2024`; no substitution of a newer library is automatic. The older available Enrichr calls did not pass a custom background, so copying those calls would not implement the final manuscript’s background statement.

For each test in the offline mode, let N be the measured background size, n the query size, M the pathway size after intersection with background, and k the overlap. Compute `phyper(k-1, M, N-M, n, lower.tail=FALSE)`. Empty queries and zero overlaps have p=1. These tests have been checked against a separate Fisher exact calculation. No Enrichr combined score is invented.

Display filtering occurs after adjustment. The available source excludes broad parent pathways and selects up to three terms per priority category by the best raw p-value across groups: common to all three; both HFpEF groups only; otherwise Control-containing; higherEF only; lowerEF only. “Control-containing” can include a term shared with one HFpEF group, so it is not called “Control only”. Dot size is overlapping gene count, colour is -log10(FDR), and grey indicates no FDR-significant result. If no pathways pass FDR, the output states that explicitly. An empty/missing term in a cached service result is not assigned an invented p-value.

### Enrichr exports

Set `enrichment_mode = "enrichr_exports"`. Place the following under the chosen data directory's `enrichr_exports/` subfolder:

- `Control.tsv`, `higherEF.tsv`, `lowerEF.tsv`: unfiltered service results with `Term`, `P.value`, `Adjusted.P.value` and `Overlap` (`k/M`) or semicolon-delimited `Genes`. Space-separated standard Enrichr header variants are accepted.
- `manifest.csv`: one row per group, with `Group`, `Library`, `BackgroundMD5`, `QueryMD5`. The library and hashes must match those exported by this workflow for the actual analysed data. The manifest documents the query and background declared for the original service run; a matching hash alone cannot prove what was submitted historically.

The export mode writes the current gene-list files before validating the manifest. Compare the actual archived submission lists; do not fill a manifest from current hashes merely to bypass a mismatch.

### Optional live Enrichr

The maintained `enrichR` interface supports a custom `background` argument; the underlying service also exposes a custom-background endpoint. See the [official Ma’ayan Lab client](https://maayanlab.github.io/maayanlab-bioinformatics/_modules/maayanlab_bioinformatics/api/speedrichr.html).

1. Use authorised study inputs with genuine gene symbols and a verified measured background.
2. Install the optional dependency: `Rscript install_dependencies.R --enrichr`.
3. Set `enrichment_mode = "enrichr_live"` in `config.R`.
4. Run `Rscript run_all.R --data YOUR_DATA_DIRECTORY --out YOUR_OUTPUT_DIRECTORY --allow-enrichr-network`.

This mode transmits the query and background **gene lists** to Enrichr. It does not send participant identifiers, medication values or expression matrices. The default synthetic workflow never uses this mode. Live service execution has not been tested. Validation covered local parsing, provenance checks and the complete offline workflow. Server-adjusted p-values are retained rather than recomputed from a filtered response. Results are cached for subsequent offline plotting.

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
