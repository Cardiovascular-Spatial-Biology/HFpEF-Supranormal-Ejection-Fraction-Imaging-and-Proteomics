# Input data dictionary

All supplied files are under `data/synthetic/`. Identifiers starting `SYN_` and pathways starting `DEMO_` are invented. No file contains study participant data.

| File | Rows / columns | Meaning |
|---|---|---|
| `expressiondata_restandstress.csv` | 480 assays × 96 samples, plus `SeqId` | Positive fictional RFU values. The first column is the unique assay key. All other headers are `SampleID` values. |
| `expressiondata-rest.csv` | 480 assays × 48 participants, plus `SeqId` | Redundant resting subset; headers are `PatientID`. Validated against the combined matrix. |
| `expressiondata-stress.csv` | 480 assays × 48 participants, plus `SeqId` | Redundant exercise subset, in the same format. |
| `sampledata.csv` | 96 samples | Explicit pairing: `SampleID`, `PatientID`, `Condition` (`Rest` or `Stress`). Suffix parsing is not used. |
| `phenotypedata.csv` | 48 participants | Minimal synthetic group, adjustment and imaging/clinical variables; described below. |
| `featuredata.csv` | 480 assays | `SeqId`, `GeneSymbol`, `TargetName`, `UniProtID`, `Role`. Multiple UniProt keys use `|`. Several assays share one gene label. |
| `proteinatlas.tsv` | 408 fictional annotation entries | `Uniprot`, `RNA.tissue.specific.nTPM`, `Reference`. The nTPM field is `tissue: value;tissue: value`. This is not a real Human Protein Atlas export. |
| `pathways.tsv` | One row per pathway–gene membership | `Pathway`, `GeneSymbol`, `Family`, `Level`, `Reference`. Thirty fictional sets at invented levels 2 or 3; not Reactome. |

`SYNTHETIC_DATA.txt` identifies the dataset and its generator seed. It is not a patient data file.

## Phenotype fields

| Field | Type / intended units | Role |
|---|---|---|
| `PatientID` | Unique string | Participant key |
| `Phenotype` | `Control`, `higherEF`, `lowerEF` | Analysis group; 16 synthetic participants per group |
| `Age` | Numeric years | Adjustment variable |
| `Sex` | `Female` / `Male` | Adjustment variable |
| `BMI` | Numeric kg/m² | Adjustment variable |
| `DeltaLVSVi` | Numeric mL/m² | Exercise minus rest LV stroke volume index |
| `DeltaLVEF` | Numeric percentage points | Exercise minus rest LVEF |
| `DeltaCI` | Numeric L/min/m² | Exercise minus rest cardiac index |
| `DeltaLAESVi` | Numeric mL/m² | Exercise minus rest LA end-systolic volume index |
| `DeltaLungWater` | Numeric percent change | Lung-water change, supplied explicitly |
| `SixMinuteWalk` | Numeric metres | Six-minute walk distance |
| `KCCQ` | Numeric score | Clinical summary score |
| `ACEi`, `ARB`, `MRA`, `BetaBlocker`, `SGLT2i`, `LoopDiuretic` | Numeric 0=no use, 1=use; NA permitted for local study inputs | Six separately modelled medication classes for Table S14 |

These are invented demonstration values with illustrative units, not estimates of the study's distributions. A few imaging/clinical values are intentionally missing to demonstrate variable-specific complete-case analysis. The expression matrices and adjustment fields have no missing values.

Medication values are independently generated using a separate random stream. They do not derive from study prescribing patterns, and no medication effect was inserted into the expression matrices. Both medication strata are represented within each synthetic group. Real missing values are not converted to zero; models use medication-specific complete cases and report their sizes.

The prefix `SYN` abbreviates **synthetic**. `P` identifies a participant, `A` an assay, `G` a gene and `U` a fictional UniProt-style key. These are labels, not clinical categories or codes for real records.

Baseline cystatin is calculated as `log2(resting RFU)` from the single assay whose feature `Role` is `cystatin_covariate`. Other features have `Role=protein`. This label only emulates the original assay-based adjustment; the synthetic target is not a real cystatin reagent.

## Supplying authorised study inputs locally

The workflow uses the explicit schema defined above. Legacy matrices with empty first-column names, space-containing clinical headers or implicit sample pairing are **not drop-in inputs**. Map them to these fields using an audited local adapter; do not silently overwrite IDs or infer a match from row order. Keep authorised study inputs outside this public package, or in ignored `data/private/`.

Use `Rscript run_all.R --data YOUR_DATA_DIRECTORY --out YOUR_OUTPUT_DIRECTORY` after matching the schema. All files are matched by IDs; row/column order may differ. Unknown or duplicated participant/sample/assay IDs fail explicitly. Positive RFU, one sample of each condition per participant, full-rank models and adequate residual degrees of freedom are checked.

For study CAMERA results, replace fictional pathways with an explicit versioned Reactome membership export and justified family assignments. The code filters the supplied `Level` field; it does not calculate hierarchy depth or download a database. Confirm the release, hierarchy convention, gene naming and measured background with the analysis authors. Clinical delta definitions and units must also be verified during conversion.

For study Fig 4d, choose `enrichr_exports` or `enrichr_live` in `config.R`. Offline fictional enrichment is deliberately restricted to `SYN_` genes and `DEMO_` pathways. Optional exports require three group TSV files and a manifest matching the current library, query lists and measured background; see `FIG4D_AND_TABLES14.md`. No participant values or expression matrix are needed by Enrichr.
