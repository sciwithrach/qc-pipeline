# utricleqc

QC pipeline for mouse utricle single-cell data: cell calling (CellBender),
doublet detection (vaeda + scDblFinder), minimal filtering, dimensionality
reduction (TopOMetry), clustering (scSHC + clustree) and outlier removal.

This folder is laid out following nf-core directory conventions
(`main.nf`, `workflows/`, `modules/local/<name>/main.nf`, `conf/*.config`,
`nextflow_schema.json`), ported from the original flat pipeline in the
parent `nextflow/` directory. `nf-core pipelines lint` has **not** been run
against it (the nf-core CLI wasn't available when this was created) — treat
the layout as lint-oriented, not lint-verified.

## Two entry points

This repo contains two pipelines that share modules, not one pipeline with
two modes, so they have two separate entry scripts:

1. **`main.nf`** — the primary QC pipeline (cell calling → doublets →
   minimal filter → TopOMetry → clustering). Run with:
   ```bash
   nextflow run main.nf -profile test
   ```
2. **`workflows/outlier_removal.nf`** — takes the `adata_doublets.h5ad`
   output of a previous `main.nf` run and removes outliers, re-running
   TopOMetry/clustering on the filtered result. 
   Run it directly:
   ```bash
   nextflow run workflows/outlier_removal.nf -profile alldata
   nextflow run workflows/outlier_removal.nf -profile postnatal
   ```

## Requirements

- [Nextflow](https://www.nextflow.io/) >= 24.10
- [Apptainer/Singularity](https://apptainer.org/) — all processes run in
  prebuilt containers (see `containers/`) which can be built as follows:
  ```bash
  bin/helpers/setup-containers.sh
  ```

## Parameters

See `nextflow_schema.json` for the full parameter list. Key ones:

| Param | Used by | Description |
|---|---|---|
| `raw_adata` | `main.nf` | Input AnnData for the QC pipeline |
| `batch` | both | obs column with batch IDs |
| `celltype_csv` | both | CSV for simple marker-based cell typing |
| `metadata` | `main.nf` | obs column used for QC plot colouring |
| `skip_cellbender` | `main.nf` | Skip CellBender (default: true) |
| `skip_scSHC` | both | Skip scSHC clustering (memory intensive) |
| `adata` | `outlier_removal.nf` | Doublet-called AnnData to filter |
| `clusters`, `doublets`, `descriptor` | `outlier_removal.nf` | Outlier-detection inputs |

## Profiles

- `test` — small bundled test data (`test/adata_test.h5ad`), for `main.nf`.
- `myriad` — UCL Myriad SGE cluster.

## Note on Cellbender

This workflow currently only runs Cellbender without GPU acceleration.

As such, I recommend you run Cellbender separately first with GPU acceleration
and use the `skip_cellbender` parameter to use your output in this pipeline.

## License

Released under the standard MIT license.
