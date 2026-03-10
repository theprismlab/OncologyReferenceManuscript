# PRISM Oncology Reference Manuscript

Data processing and biomarker analysis pipeline for the PRISM Oncology Reference dataset.

## Scripts

| Script | Description |
|--------|-------------|
| `UTILITIES.R` | Shared functions: dose-response curve fitting, univariate biomarker analysis, target recovery, and random forest biomarker models |
| `1 - DATA_PROCESSING.R` | Data processing pipeline: noise floor correction, spline normalization, QC, log-fold-change computation, replicate collapsing, dose-response curve fitting, and output matrix generation |
| `2 - DATA_PROCESSING - UNCORRECTED FILES FOR COMPARISON.R` | Same as script 1 but skips the cell-line artifact regression step, producing uncorrected outputs for benchmarking |
| `3 - GENERATE BIOMARKER TABLES.R` | Computes univariate biomarkers, fits random forest predictive models with cross-validation, and generates summary scores (polypharmacology, selectivity, variable importance) |

## Required Input Data

Scripts expect the following directory structure relative to the project root:

```
data/
├── input data/
│   ├── PRISMOncologyReferenceInstMeta.csv
│   ├── PRISMOncologyReferenceAnalyteMeta.csv
│   ├── PRISMOncologyReferenceLMFI.csv
│   └── PRISMOncologyReferenceCompoundList.csv
└── external inputs/
    ├── depmap_oncref_manuscript.h5
    └── TKRTK_gene_symbols.csv
```

Outputs are written to `data/processed data/` and `results/`.

## Running with Docker

### Build the image

```bash
./docker_build.sh
```

### Run the full pipeline

Pass the path to your data directory as an argument. The script creates output subdirectories and runs all three scripts in sequence:

```bash
./docker_run.sh /path/to/data
```

The data directory should contain the input files described above. Processed outputs and results will be written into subdirectories within the same data directory.

### Resource considerations

Scripts 1 and 2 use `parallel::mclapply` with `detectCores() - 1` for dose-response curve fitting. The biomarker pipeline (script 3) is memory-intensive due to large matrix operations. You can control CPU and memory allocation:

```bash
./docker_run.sh /path/to/data --cpus 8 --memory 16g
```
