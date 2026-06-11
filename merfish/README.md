# MERFISH CLI — batch spatial-transcriptomics pipeline

A SLURM-ready, text-config-driven command-line port of the PhenoSuite
**MERFISH** module (`phenosuite/merfish/app.R`). It runs the full
segmentation-based spatial-transcriptomics workflow headlessly, one sample per
SLURM **array task**, so an arbitrary number of samples process in parallel from
a single submission.

```
Import → QC → Normalize → HVG → PCA → UMAP → Leiden cluster →
Neighborhood enrichment → Spatially variable genes → Cluster markers → Export
```

## Layout

| File | Role |
|------|------|
| `merfish-config.txt`   | **Edit this.** All paths, column mappings, and parameters. |
| `batch-inputs/`        | Line-aligned per-sample manifests (row *N* = sample *N*). |
| `merfish-meta.s`       | Entry point — validates manifests, submits the array job. |
| `run-merfish.s`        | SLURM array task — extracts row *N*, calls the R driver. |
| `run-merfish.R`        | CLI entry point (arg parsing, validation). Standalone-runnable. |
| `RunMerfish.R`         | Pipeline orchestrator (`RunMerfish()`). |
| `merfish-utils.R`      | Analysis primitives (QC, normalize, PCA, UMAP, Leiden, spatial stats, DE, figures). |
| `environment.yml`      | Conda env (`runmerfish`). |
| `makeRunLog-batch.sh`  | Snapshot config + manifests to `run-logs-batch/`. |

## The manifests (one line per sample)

Each list in `batch-inputs/` holds **one row per sample**; row *N* across every
list is the same sample. Pad any missing per-sample value with the literal
string `NULL`.

- `expression_files.txt` — cell × gene matrix per sample (CSV/TSV/`.gz`; first column = cell ID)
- `metadata_files.txt` — per-cell metadata per sample (first column = cell ID; must contain the X/Y coordinate columns named in the config)
- `out_dirs.txt` — output directory per sample
- `sample_ids.txt` — output filename prefix per sample (optional; `NULL` → output dir basename)

`batch_size` is auto-computed from the expression list (`wc -l`), so the array is
sized automatically.

## Quick start

```bash
# 1. One-time: provision the R environment
conda env create -f environment.yml
conda activate runmerfish

# 2. Point the manifests at your data and tune parameters
$EDITOR merfish-config.txt
$EDITOR batch-inputs/expression_files.txt batch-inputs/metadata_files.txt batch-inputs/out_dirs.txt

# 3. (optional) Archive this submission's config for reproducibility
bash makeRunLog-batch.sh

# 4. Submit the whole batch as a SLURM array
sbatch merfish-meta.s
```

`merfish-meta.s` checks that every manifest has at least `batch_size` rows, then
submits `run-merfish.s` as `--array=1-<batch_size>`. Each task writes its SLURM
logs to `RunMerfish_<jobid>_<taskid>.{out,err}`.

## Run one sample directly (no SLURM)

`run-merfish.R` is fully standalone — handy for testing on a login/dev node:

```bash
conda activate runmerfish
Rscript run-merfish.R \
  --expression-file=/path/sampleA/cell_by_gene.csv.gz \
  --metadata-file=/path/sampleA/cell_metadata.csv.gz \
  --out-dir=/path/results/sampleA \
  --sample-id=sampleA \
  --x-col=center_x --y-col=center_y \
  --area-col=cell_area --negctrl-col=blank_counts

Rscript run-merfish.R --help    # full option list
```

## Outputs (per sample, prefixed with `<sample_id>_`)

| File | Contents |
|------|----------|
| `metadata_clusters.csv`       | Cell metadata + `cluster`, `Phenotype`, UMAP coords |
| `normalized_expression.csv`   | Normalized cell × gene matrix |
| `cluster_markers.csv`         | One-vs-rest Wilcoxon markers per cluster |
| `svg.csv`                     | Spatially variable genes (Moran's I, BH-adjusted p) |
| `neighborhood_enrichment.csv` | Cluster × cluster co-localization Z-scores |
| `spe.rds`                     | `SpatialExperiment` (assays `exprs`/`logcounts`/`counts`, `UMAP`/`PCA` reducedDims) for downstream PhenoSuite modules |
| `figures/*.pdf`               | QC violin, spatial clusters, UMAP clusters |
| `run_manifest.json`           | Inputs, parameters, counts, timing, git SHA — full provenance |

The `spe.rds` carries a `Phenotype` column so it drops straight into the
downstream phenotyping / PCF / spatial-interaction modules.

## Notes

- **Expression orientation:** input is read as cells × genes. If your matrices
  are genes × cells, set `transpose=TRUE` in the config (or pass `--transpose`).
- **Disabling a QC filter:** set its config value to `NULL`. Area and density
  filters additionally require `area_col`; the negative-control filter requires
  `negctrl_col`.
- **Tuning SLURM:** `module1_cpus`, `module1_mem`, `module1_time`, and
  `module1_partition` in the config control the per-task resource request.
- **`module load r/4.1.2`** in `run-merfish.s` is site-specific; on a system
  using the conda env instead, remove it and `conda activate runmerfish` before
  `sbatch`, or add the activation into `run-merfish.s`.
