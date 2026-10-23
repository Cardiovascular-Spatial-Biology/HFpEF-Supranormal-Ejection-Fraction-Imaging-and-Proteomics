# Validation

The synthetic workflow completed on 2026-10-23, producing 11 figures/tables. All 16 automated checks passed.

Checks cover reproducible data generation, ID matching, paired-model estimates, duplicate-ID rejection, assay mapping, PLS-DA calculations, CAMERA FDR and medication contrasts. These checks validate the demo workflow, not reproduction of the manuscript results.

Run from the project folder:

```sh
Rscript run_all.R
Rscript tests/validate_demo.R
```

The CAMERA output check reads `results/demo/`. Tested versions are in [software_versions.csv](software_versions.csv); distributed-file hashes are in [SHA256SUMS.txt](SHA256SUMS.txt).
