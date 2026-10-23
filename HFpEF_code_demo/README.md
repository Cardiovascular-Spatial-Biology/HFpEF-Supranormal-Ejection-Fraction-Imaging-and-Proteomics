# HFpEF imaging and proteomics

Analysis code for **HFpEF with supranormal ejection fraction is a distinct exercise endotype and a candidate for cardiac myosin inhibition**.

This R workflow includes synthetic data for running the proteomics analyses. All supplied participant records are fictional. 

Some modules preserve verified analytical code. Others were reconstructed from the available scripts and written methods and have not been matched to the final manuscript results. Their status is listed below and in [Methods and provenance](docs/METHODS_AND_PROVENANCE.md).

## Quick start

In RStudio, open `HFpEF_code_demo.Rproj`, then run:

```r
# Install missing dependencies once; internet access is required.
source("install_dependencies.R")

# Run the analysis in a fresh R session.
source("run_all.R")
```

The default run reads `data/synthetic/` and writes figures and tables to `results/demo/`. It runs offline after dependencies are installed. If that output folder is already populated, a subsequent run creates a timestamped folder under `results/`.

Alternatively, run these commands from the project folder:

```sh
Rscript install_dependencies.R
Rscript run_all.R
Rscript tests/validate_demo.R
```

The validation script checks the supplied inputs and analytical functions. Its CAMERA output check uses `results/demo/`, so run the workflow first. See [Validation](docs/VALIDATION.md) for the checks performed.

To choose a different output folder:

```sh
Rscript run_all.R --out results/my_run
```

An explicitly selected output folder must be empty or new. Input files are not overwritten, and no personal working-directory setting is needed.

## Included analyses

| Module | Analysis | Figure or table; status |
|---|---|---|
| `R/01_paired_limma.R` | Paired exercise response within each group, adjusted for age, sex, BMI and baseline cystatin | Fig 4b analytical core; verified against the available study results |
| `R/02_overlap_ternary.R` | Overlap counts, Venn diagram and ternary plot | Fig 4b-c; regenerated layout |
| `R/03_tissue_plots.R` | Tissue-associated exercise responses with corrected assay-ID mapping | Fig 4e-f; regenerated layout |
| `R/04_resting_comparisons.R` | Unadjusted resting HFpEF versus control and higherEF versus lowerEF comparisons | Fig S4; reconstructed |
| `R/05_pca_plsda.R` | Resting PCA and two PLS-DA displays, with corrected X-variance labels | Fig S5a-c |
| `R/06_camera_pathways.R` | Associations between protein changes and imaging/clinical variables | Fig 4g; reconstructed; final study reference files are not supplied |
| `R/07_pathway_enrichment.R` | Measured-background enrichment and grouped dot plot | Fig 4d; reconstructed; fictional pathways in the default demonstration |
| `R/08_medication_sensitivity.R` | Six medication classes, each with a main effect and an exercise interaction | Table S14; reconstructed; the exercise contrast requires confirmation against the final author model |

`R/00_common.R` provides input checks, ID matching and shared functions. Analysis settings are in `config.R`; `run_all.R` loads and runs the modules.

Fig 4d uses an offline hypergeometric test on fictional pathways by default. Optional modes read Enrichr exports or submit gene lists to Enrichr with a measured-gene background. Table S14 reports medication-adjusted average exercise responses and interaction tests. Definitions, assumptions and optional modes are documented in [Fig 4d and Table S14](docs/FIG4D_AND_TABLES14.md).

Clinical-imaging statistics, pressure-volume modelling and digital-twin simulations are outside the scope of this package. This folder is not a complete reproduction of every analysis in the manuscript.

## Data

The example contains **48 fictional participants, 96 paired samples, 480 assays and 400 gene labels**. Eight CSV/TSV input files provide expression matrices, sample metadata, phenotype fields, feature annotations and reference sets. File formats and units are described in the [Data dictionary](docs/DATA_DICTIONARY.md).

`SYN` stands for **synthetic**. For example, `SYN_P001` identifies fictional participant 001. These identifiers have no link to real participants. Medication fields are also synthetic.

To recreate the supplied data in a new folder:

```sh
Rscript generate_demo.R --out data/synthetic_regenerated
```

The generator uses seed `20261019` and reads no study files. It refuses to overwrite a nonempty destination. Neither the expression values nor the reference annotations were derived by modifying patient records.

## Outputs and reproducibility

Figures are saved as PDF and PNG; model results, mappings and sample counts are saved as CSV. Demonstration figures are labelled as synthetic. Each run records input checksums and software versions without personal paths or a full session dump.

[Software versions](docs/software_versions.csv) records the tested environment; it is not a dependency lockfile. [SHA256SUMS.txt](docs/SHA256SUMS.txt) records the distributed files. [distribution_check.json](docs/distribution_check.json) summarises the file-content review.

The distribution contains code, synthetic inputs, documentation and tests. Generated results, private data directories, local configuration and session caches are excluded by `.gitignore`. Review any additional files before publishing them; ignored files can still be added explicitly to Git.
