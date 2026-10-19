# Validation record

Validated on 2026-10-19. These checks cover the synthetic demonstration and the analytical components described below; they do not establish reproduction of every manuscript analysis.

- The complete eight-step workflow, including Fig 4d and Table S14, ran in a fresh `Rscript` process on the supplied synthetic inputs and produced all implemented figures and CSV outputs.
- The complete workflow also ran when launched from an unrelated working directory with an absolute script location and an output directory containing spaces. The package has no personal working-directory settings.
- The generator reproduced every supplied synthetic input byte using its fixed seed.
- Input rows and expression columns were independently shuffled. Explicit ID matching preserved paired values and the group-wise limma estimates.
- Duplicate participant IDs were rejected explicitly.
- All 488 expanded fictional UniProt mappings retained valid assay IDs.
- CAMERA selected exactly one of 480 assays for each of 400 gene labels according to highest mean resting RFU.
- PLS-DA axis percentages agreed with an independent projection calculation from the scaled X matrix and scores.
- CAMERA BH-adjusted values agreed with a separate per-variable `p.adjust` calculation.
- The reusable validation script passed 23 checks. Run `Rscript tests/validate_demo.R` after a full demonstration run to include the CAMERA output check.
- Offline Fig 4d enrichment agreed with an independent Fisher exact calculation. Tests covered empty queries, zero-overlap sets, background intersection, unique-gene counting and preserving supplied Enrichr adjusted p-values. Live Enrichr was blocked without its explicit flag; no live service analysis was performed.
- The medication baseline reproduced the original paired model. The centred exercise coefficient and moderated t agreed with an independently constructed prevalence-weighted contrast from an uncentred model. Reversing exposure coding preserved the average response and reversed the interaction sign. Constant exposure was correctly marked not estimable. Interaction BH correction used all tested assays.
- In a separate private validation, the paired limma function was run on the authorised study inputs and compared with previously reproduced group-wise results. Maximum absolute differences across logFC, mean expression, moderated t, raw p, adjusted p and B statistics were below `6e-14`. No study matrix, participant record or private adapter is included in this package.
- Representative generated figures were visually inspected for labels, legends and readable output. Layout is regenerated and is not asserted to be pixel-identical to the manuscript.

The paired-model equivalence check does not validate reconstructed S4, CAMERA, Fig 4d or Table S14 against the final author's outputs. Their limitations and explicit implementation choices are documented in `METHODS_AND_PROVENANCE.md` and `FIG4D_AND_TABLES14.md`.

A separate content review found no Chinese characters in filenames or text, personal paths, known study identifier patterns, email addresses or credential patterns. All supplied inputs were reproduced byte-for-byte by the generator. The distribution contains 35 text files: code, synthetic inputs, documentation, tests and an empty results placeholder. It contains no Git history, cached R objects, logs or generated results. These findings apply to this folder and archive; the existing online repository and its history were not examined.
