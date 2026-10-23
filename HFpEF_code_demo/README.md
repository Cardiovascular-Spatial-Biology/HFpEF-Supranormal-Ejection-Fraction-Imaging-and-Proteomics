# HFpEF imaging and proteomics

The supplied data are synthetic: 48 fictional participants, 96 paired samples and 480 assays. Demo outputs are not the study results. See [Methods](docs/METHODS_AND_PROVENANCE.md) for implementation details and limitations.

## Run

From this folder:

```sh
Rscript install_dependencies.R
Rscript run_all.R
Rscript tests/validate_demo.R
```

In RStudio, open `HFpEF_code_demo.Rproj`, then use `source("install_dependencies.R")` and `source("run_all.R")` in a fresh session. Installation requires internet access; the analysis runs offline.

Inputs are in `data/synthetic/`. Results are saved to `results/demo/` as PDF, PNG and CSV files. Later runs use a timestamped folder if `results/demo/` is populated. The validation command checks CAMERA outputs in `results/demo/`.

## Scripts

| Script | Analysis | Figure/table |
|---|---|---|
| `R/01_paired_limma.R` | Paired exercise-response models | Fig 4b core analysis |
| `R/02_overlap_ternary.R` | Venn diagram and ternary plot | Fig 4b-c |
| `R/03_tissue_plots.R` | Tissue-associated responses | Fig 4e-f |
| `R/04_resting_comparisons.R` | Resting group comparisons | Fig S4 |
| `R/05_pca_plsda.R` | PCA and PLS-DA | Fig S5a-c |
| `R/06_camera_pathways.R` | Protein-change associations with imaging/clinical variables | Fig 4g |
| `R/08_medication_sensitivity.R` | Medication sensitivity models | Table S14 |

Settings are in `config.R`; shared input checks are in `R/00_common.R`.

## Data and options

`SYN` means synthetic. All supplied expression values, participant records and reference annotations are fictional. File formats are listed in the [Data dictionary](docs/DATA_DICTIONARY.md).

```sh
Rscript run_all.R --data YOUR_DATA_DIRECTORY --out YOUR_OUTPUT_DIRECTORY
Rscript generate_demo.R --out data/synthetic_regenerated
```

Custom inputs must match the documented schema. Output and regenerated-data folders must be empty or new. Keep real participant data outside this public package.

[Validation](docs/VALIDATION.md) · [Software versions](docs/software_versions.csv) · [File checksums](docs/SHA256SUMS.txt)
