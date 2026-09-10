# qc-pipeline

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.0-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)

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

## Running with `skip_cellbender = true` (default)

When `skip_cellbender` is `true`, `main.nf` skips cell calling entirely —
the Cellbender/batch-splitting/joining processes never run, and `raw_adata`
is used as-is as the starting point for QC. This means:

- **`raw_adata`** must already be a *cell-called* AnnData (empty droplets /
  ambient RNA already removed), not the fully unfiltered matrix — e.g. the
  output of a Cellbender run done separately with GPU acceleration (see
  above), or already filtered by your own thresholding. No empty-droplet
  removal happens implicitly in this mode.
- **`batch`** must name an existing column in `raw_adata.obs` (e.g. `"plate"`
  in `conf/test.config`) — used to group cells for batch-aware QC
  plots/reports and passed to `run_vaeda`/`run_scDblFinder` for batch-aware
  doublet detection.
- **`metadata`** must name an existing column in `raw_adata.obs` (e.g.
  `"sample"` in `conf/test.config`) — used as the colouring/grouping
  variable in QC plots and summary stats. Any numeric or categorical obs
  column works.
- **`celltype_csv`**, if provided, must be a CSV with **no header row**:
  column 1 is a cell-type label, the remaining columns are marker gene
  symbols for that type (rows can have different numbers of marker
  columns — see `test/celltype_test.csv` for a real example). Marker genes
  not found in `adata.var_names` are silently ignored, so check gene
  symbols match your `.var_names` naming convention.

If `skip_cellbender = false` instead, `raw_adata` should be the *fully
unfiltered* matrix, and the pipeline additionally accepts
`params/cellbender_params.csv` — an optional per-batch settings file
(columns: `batch, total_droplets_included, expected_cells`; omit a value or
the whole file to let Cellbender estimate automatically).

## License

Released under the standard MIT license.
