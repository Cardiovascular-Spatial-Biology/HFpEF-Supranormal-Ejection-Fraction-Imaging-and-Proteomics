# Input data dictionary

All files below are in `data/synthetic/` and contain synthetic values. `SYN_` identifiers and `DEMO_` pathways are fictional.

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

`SYNTHETIC_DATA.txt` records the generator seed.

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

Units are illustrative. Some imaging/clinical values are missing; RFU and adjustment fields are complete. Medication values are independently generated with no medication effect inserted into RFU.

`SYN_P`, `SYN_A`, `SYN_G` and `SYN_U` identify fictional participants, assays, genes and UniProt-style keys. Baseline cystatin is `log2(resting RFU)` for the assay with `Role=cystatin_covariate`.

## Custom inputs

Use the same field names and explicit IDs; sample and assay order may differ. Provide positive RFU and one Rest and one Stress sample per participant. Medication NA values are handled as missing, not zero.

```sh
Rscript run_all.R --data YOUR_DATA_DIRECTORY --out YOUR_OUTPUT_DIRECTORY
```

Keep real participant data outside the public package. For CAMERA, provide versioned pathway membership, hierarchy levels and family assignments. The code uses the supplied `Level` field and does not download reference data. Verify clinical units and delta definitions before use.
