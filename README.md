# phenosuite-CLI

Command-line interface for the core PhenoSuite pipelines, designed to run on SLURM-based HPC clusters. The repository bundles three independent spatial-biology modules that operate on multiplexed tissue imaging data (CODEX, multiplexed IF, HALO / Mesmer segmentations):

- **[RunPhenomenalist](RunPhenomenalist/)** — cellular scaling, dimensionality reduction, and clustering from single-cell segmentation tables (R).
- **[masquerade](masquerade/)** — circular cluster-mask generation for QuPath overlays on OME-TIFF / QPTIFF images (Python).
- **[spatial-dynamics](spatial-dynamics/)** — pairwise log-odds and multi-cell-type neighborhood enrichment for spatial cell–cell relationships (Python).

Each module is a self-contained `sbatch`-driven array job: a config file declares parameters, line-aligned text files in `batch-inputs/` declare the samples, and one `sbatch …` command fans the batch out across array tasks.

---

## Table of contents

1. [Repository layout](#repository-layout)
2. [Prerequisites](#prerequisites)
3. [Quickstart](#quickstart)
4. [Modules](#modules)
   - [RunPhenomenalist](#runphenomenalist)
   - [masquerade](#masquerade)
   - [spatial-dynamics](#spatial-dynamics)
5. [Configuration reference](#configuration-reference)
6. [Active vs. legacy files](#active-vs-legacy-files)
7. [Run logs](#run-logs)
8. [Known limitations](#known-limitations)
9. [Troubleshooting](#troubleshooting)

---

## Repository layout

```
phenosuite-CLI/
├── README.md
├── .gitignore
│
├── RunPhenomenalist/                         # R pipeline: scaling + dimensionality reduction + clustering
│   ├── RunPhenomenalist.R                    #   core pipeline function
│   ├── run-phenomenalist.R                   #   CLI wrapper (parses 9 positional args)
│   ├── phenomenalist-utils.R                 #   plotting / clustering / cell-type utilities
│   ├── phenomenalist-config.txt              #   batch config (edit this)
│   ├── phenomenalist-meta.s                  #   SLURM entry point — `sbatch` this
│   ├── run-phenomenalist.s                   #   SLURM array task (called by meta)
│   ├── makeRunLog-batch.sh                   #   captures config + input lists per run
│   ├── batch-inputs/                         #   line-aligned per-sample inputs
│   │   ├── segmentation_files.txt
│   │   ├── out_dirs.txt
│   │   ├── failed-markers.txt
│   │   ├── nuclear-markers.txt
│   │   └── labels.txt
│   ├── run-logs-batch/                       #   timestamped run logs (generated)
│   └── v0/                                   #   legacy docopt interface (archived)
│
├── masquerade/                               # Python pipeline: cluster-mask OME-TIFFs
│   ├── Masquerade.py                         #   core class (TIFF I/O, mask geometry)
│   ├── masquerade_interface.py               #   CLI wrapper (11 positional args)
│   ├── masquerade_utils.py                   #   legacy standalone helpers (superseded)
│   ├── Masquerade_v0.py                      #   legacy function-based impl (archived)
│   ├── masquerade_interface_v0.py            #   legacy wrapper (archived)
│   ├── masquerade_interface_v1.py            #   experimental wrapper (not used by launcher)
│   ├── environment.yml                       #   conda env spec (python=3.10 + deps)
│   ├── configFile-batch.txt                  #   batch config (edit this)
│   ├── run-masquerade.s                      #   SLURM entry point — `sbatch` this
│   ├── run-masquerade-batch.sh               #   SLURM array task (calls python + bfconvert)
│   ├── run-masquerade-batch_v0.sh            #   legacy launcher (archived)
│   ├── makeRunLog-batch.sh                   #   captures config + input lists per run
│   ├── bftools/                              #   bundled Bio-Formats CLI (bfconvert)
│   ├── batch-inputs/                         #   line-aligned per-sample inputs
│   │   ├── qptiff-batch.txt
│   │   ├── spatial_anno-batch.txt
│   │   ├── outPaths-batch.txt
│   │   └── marker-metadata-batch.txt
│   └── run-logs-batch/                       #   timestamped run logs (generated)
│
└── spatial-dynamics/                         # Python pipeline: PWLO + neighborhoods
    ├── pwlo_es_pt.py                         #   pairwise log-odds implementation
    ├── run-pwlo.py                           #   CLI wrapper for PWLO
    ├── n_simplex_neighborhoods.py            #   n-way (3+ celltype) neighborhood enrichment
    ├── getNeighborhoods.py                   #   wrapper orchestrating neighborhood analysis
    ├── run-spatial_circuit-enrichment.py     #   CLI wrapper for circuit enrichment
    ├── optimize-neighborhoods.py             #   threshold optimization + animation
    ├── config-spatial_dynamics.txt           #   batch config (edit this)
    ├── run-spatial_dynamics.s                #   SLURM entry point — `sbatch` this
    ├── run-pwlo.s                            #   SLURM array task (module=0)
    ├── run-spatial_circuit-enrichment.s      #   SLURM array task (module=1)
    ├── batch-inputs/                         #   line-aligned per-sample inputs
    │   ├── spatial-annos.txt
    │   ├── outs.txt
    │   └── labels.txt
    └── *-v0.R, *-v1.R, optimize-neighborhoods.R   # legacy R implementations (archived)
```

---

## Prerequisites

Only the **masquerade** module ships a dependency manifest ([`masquerade/environment.yml`](masquerade/environment.yml)). The R and spatial-dynamics environments must be provisioned by hand using the lists below.

### Cluster

- A **SLURM** cluster. Partition names in the shipped configs (`a100_short`, `cpu_dev`) are site-specific — edit them for your cluster before running.

### RunPhenomenalist

Provision the shipped conda env and install the public [phenomenalist](https://github.com/igordot/phenomenalist) package from GitHub:

```bash
cd RunPhenomenalist/
conda env create -f environment.yml
conda activate runphenomenalist
Rscript -e 'remotes::install_github("igordot/phenomenalist")'
```

What [environment.yml](RunPhenomenalist/environment.yml) pins:

- `r-base=4.1` + `r-remotes`
- CRAN deps of phenomenalist: `cowplot`, `data.table`, `dplyr`, `FNN`, `ggplot2`, `ggsci`, `glue`, `igraph`, `janitor`, `RColorBrewer`, `readr`, `rlang`, `scattermore`, `stringr`, `tibble`, `tidyr`, `tidyselect`, `tidyverse`, `uwot`
- Extra CRAN deps used by `RunPhenomenalist.R` / `phenomenalist-utils.R`: `circlize`, `mclust`
- Bioconductor deps: `ComplexHeatmap`, `MatrixGenerics`, `scran`, `scuttle`, `SingleCellExperiment`, `SpatialExperiment`, `SummarizedExperiment`

The local `RunPhenomenalist/phenomenalist-utils.R` file is still sourced at runtime — it provides `create_object.mod`, `cluster.mod`, `plot_heatmap.mod_v1`, `plot_dr.mod`, `plot_spatial.mod`, `prepare_mask_inputs`, `phenomenalist.preprocess`, and `assign_celltype_with_template`, which are local modifications / helpers not in the public package.

If you're not using conda, install the R packages however you normally would and make sure `library(phenomenalist)` works before running the pipeline.

### masquerade

Create the conda env from the shipped spec — the launcher expects it to be named `masquerade` and runs `source activate masquerade`:

```bash
cd masquerade/
conda env create -f environment.yml
# or, if the env already exists:
conda env update -f environment.yml --prune
```

The spec pins `python=3.10` and pulls `tifffile`, `imagecodecs`, `scipy`, `pandas`, `matplotlib`, `scikit-image`, `tqdm`, `numpy` from `bioconda` / `conda-forge`. Note that [environment.yml](masquerade/environment.yml) contains a `prefix:` line pointing at the original author's `$HOME` — conda ignores it when the file is passed to `env create`, but you can strip it if you like a clean file.

Also required at runtime:

- **Java 17** for the bundled Bio-Formats `bfconvert` (the launcher runs `module load java/17.0.0`). `bfconvert` ships in [masquerade/bftools/](masquerade/bftools/) and is used to convert outputs ≥ 4 GB to pyramidal OME-TIFF.

### spatial-dynamics

- A Python 3 env exposed via `module load condaenvs/gpu/machinelearning` (edit to match your cluster). Required packages:
  - `numpy`, `pandas`
  - `scipy`, `scikit-learn`
  - `matplotlib`, `seaborn`

---

## Quickstart

All three modules follow the same pattern:

1. Edit the config file (`*-config*.txt`) in the module directory.
2. Add one line per sample to each file in `batch-inputs/` — **line N across every file = sample N**. The array size is auto-computed from the line count of the primary input list.
3. Submit the meta launcher with `sbatch`.

### RunPhenomenalist

```bash
cd RunPhenomenalist/
# 1. Edit phenomenalist-config.txt (HALO, clustering_res, SLURM params, etc.)
vi phenomenalist-config.txt
# 2. Populate batch-inputs/ — one line per sample in each file; use 'NULL' for empty rows
vi batch-inputs/segmentation_files.txt
vi batch-inputs/out_dirs.txt
vi batch-inputs/failed-markers.txt
vi batch-inputs/nuclear-markers.txt
vi batch-inputs/labels.txt
# 3. Submit. phenomenalist-meta.s pre-validates that every batch-inputs file
#    has >= batch_size rows before calling sbatch.
sbatch phenomenalist-meta.s
```

To run the underlying CLI on a single sample outside of SLURM (useful for debugging):

```bash
Rscript run-phenomenalist.R --help
Rscript run-phenomenalist.R \
    --segmentation-file=/data/sample1.csv \
    --out-dir=/data/sample1/out \
    --halo=T \
    --clustering-res=5:7 \
    --failed-markers=AF-Ak154,DAPI
```

### masquerade

```bash
cd masquerade/
# 1. Edit configFile-batch.txt (radius, filled, target_size, SLURM params)
vi configFile-batch.txt
# 2. Populate batch-inputs/ — one line per sample in each file
vi batch-inputs/qptiff-batch.txt
vi batch-inputs/spatial_anno-batch.txt
vi batch-inputs/outPaths-batch.txt
vi batch-inputs/marker-metadata-batch.txt
# 3. Submit
sbatch run-masquerade.s
```

### spatial-dynamics

```bash
cd spatial-dynamics/
# 1. Edit config-spatial_dynamics.txt — pick module=0 (PWLO) or module=1 (circuit)
vi config-spatial_dynamics.txt
# 2. Populate batch-inputs/ — one line per sample in each file
vi batch-inputs/spatial-annos.txt
vi batch-inputs/outs.txt
vi batch-inputs/labels.txt
# 3. Submit
sbatch run-spatial_dynamics.s
```

---

## Modules

### RunPhenomenalist

#### What it does

Reads single-cell segmentation measurements (one row per cell, one column per marker), builds a `SpatialExperiment` object, performs Leiden clustering at one or more resolutions, and emits marker-expression heatmaps, UMAP / spatial expression plots, and cluster assignments. Optionally assigns cell types using a manual gating template. The output `mask-inputs/` directory is formatted to feed directly into the **masquerade** module.

Supports two segmentation formats via the `HALO` toggle: HALO (`HALO=T`) and Mesmer (`HALO=F`).

#### Inputs

- **Segmentation files** (CSV / TSV / TSV.gz): one row per cell, columns for `x`, `y`, and per-marker intensities. Format depends on the `HALO` flag.
- Optional **phenotyping template** CSV with manual gating rules (`phenotyping_template` key in the config).

#### Outputs

Written under each `${out_dir}/out-phenomenalist/`:

- `spe.rds` — serialized `SpatialExperiment` with expression, coordinates, and cluster assignments
- `*-heatmap.png` — per-cluster marker expression heatmaps
- `UMAP-expression/` — UMAP plots colored by each marker
- `spatial-expression/` — tissue-space plots colored by each marker
- `mask-inputs/` — coordinate CSVs ready for masquerade
- `phenotyping_template_results/` — annotated cell-type plots (when a template is provided)

#### Entry point and orchestration

`sbatch phenomenalist-meta.s` reads [phenomenalist-config.txt](RunPhenomenalist/phenomenalist-config.txt), validates that every `batch-inputs/*.txt` file has at least `batch_size` rows (submission aborts with a clear error otherwise), and submits [run-phenomenalist.s](RunPhenomenalist/run-phenomenalist.s) as an array job sized `1..batch_size`. Each array task extracts its row from every `batch-inputs/*.txt` file with `sed -n "${SLURM_ARRAY_TASK_ID}p"` and invokes the CLI with named flags:

```bash
Rscript run-phenomenalist.R \
    --segmentation-file="${segmentation_file_tmp}" \
    --failed-markers="${failed_markers_tmp}" \
    --nuclear-markers="${nuclear_markers_tmp}" \
    --halo="${HALO}" \
    --out-dir="${out_dir_tmp}" \
    --clustering-res="${clustering_res}" \
    --classifier-label="${classifier_label_tmp}" \
    --max-cells="${max_cells}" \
    --phenotyping-template="${phenotyping_template}"
```

Full option reference:

```
Required:
  --segmentation-file=PATH    per-cell segmentation CSV / TSV / TSV.gz
  --out-dir=DIR               output directory (created if missing)

Optional:
  --halo=T|F                  HALO (T) vs Mesmer (F) format                 [default: T]
  --failed-markers=LIST       comma-separated markers to drop               [default: none]
  --nuclear-markers=LIST      comma-separated nuclear markers               [default: none]
  --classifier-label=LIST     comma-separated classifier labels             [default: none]
  --clustering-res=SPEC       '1,2,3' explicit list, or '5:7' range         [default: 1,2]
  --max-cells=N               cell subsample threshold                      [default: 100000]
  --phenotyping-template=PATH manual-gating template CSV                    [default: none]
  -h, --help                  show help and exit
```

Sentinels that all mean *not set* for any option value: `NULL`, `null`, `NA`, `none`, `None`, `''` (empty), `0`. That makes every row of every `batch-inputs/*.txt` file safe to pad with `NULL`.

The legacy 9-positional invocation (`Rscript run-phenomenalist.R <seg> <failed> <nuclear> <HALO> <out> <res> <label> <max> <template>`) still works for backward compatibility with any external caller, but the named-flag form is preferred.

The public [`phenomenalist`](https://github.com/igordot/phenomenalist) R package is loaded via `library(phenomenalist)` — no local sourcing of package `.R` files. If it's not installed the wrapper aborts with a clear install hint.

One environment variable controls local sibling-file discovery:

| Variable | Meaning |
|---|---|
| `PHENOMENALIST_DIR` | Directory holding `RunPhenomenalist.R` and `phenomenalist-utils.R`. Auto-detected from the script location; override only if you've split the files. |

`phenomenalist-meta.s` and `run-phenomenalist.s` both set `PHENOMENALIST_DIR` to their own directory automatically.

#### batch-inputs/ format

Every file in [RunPhenomenalist/batch-inputs/](RunPhenomenalist/batch-inputs/) holds **one row per sample**, and row N must describe the same sample across all files. Use the literal string `NULL` for a row with no value.

| File | Contents | Example row |
|---|---|---|
| `segmentation_files.txt` | absolute path to cell-measurements CSV / TSV | `/data/sample1.unmixed.qptiff_measurements.tsv.gz` |
| `out_dirs.txt` | output directory (created if missing) | `/data/sample1/out-phenomenalist/` |
| `failed-markers.txt` | comma-delimited markers to skip, or `NULL` | `AF-Ak154,DAPI` |
| `nuclear-markers.txt` | comma-delimited nuclear markers, or `NULL` | `DAPI` |
| `labels.txt` | classifier label list, or `NULL` | `PAX5` |

**To add a sample:** append one line to each file above. The config's `batch_size=$(wc -l …)` will pick up the new count automatically.

---

### masquerade

#### What it does

Reads a multi-channel tissue image (OME-TIFF or QPTIFF, routinely 4+ GB), draws a circular mask (filled disk or ring) around every cell centroid in a spatial-annotation CSV, and emits a multi-channel OME-TIFF where each channel holds the mask for one cluster. The result drops straight into QuPath via its Bio-Formats reader for overlay visualization. Large outputs are automatically converted to tiled pyramidal OME-TIFF via the bundled `bfconvert`.

#### Inputs

- **Image** — OME-TIFF or QPTIFF. The module auto-detects the format and extracts channel names from OME-XML metadata or the `Biomarker` tag.
- **Spatial annotation CSV** — must contain `x`, `y`, and `cluster` columns. One row per cell.
- **Marker metadata CSV** (optional) — if provided, its `x` column lists markers to retain in the output image; pass `None` to keep all.

#### Outputs

- `<basename>.ome.tiff` — multi-channel OME-TIFF with one mask channel per cluster, hand-written OME-XML set to `SamplesPerPixel=1` so QuPath's Bio-Formats reader opens it correctly.
- `<basename>.pyramidal.ome.tiff` — only for outputs ≥ 4 GB; generated by `bfconvert` with `-tilex 512 -tiley 512 -pyramid-resolutions 5 -pyramid-scale 2 -compression LZW`. The intermediate non-pyramidal file is removed on success (see [run-masquerade-batch.sh:44-67](masquerade/run-masquerade-batch.sh)).

#### Entry point and orchestration

`sbatch run-masquerade.s` reads [configFile-batch.txt](masquerade/configFile-batch.txt) and submits [run-masquerade-batch.sh](masquerade/run-masquerade-batch.sh) as an array job sized to the line count of the spatial-metadata list. Each task activates the `masquerade` conda env, loads `java/17.0.0`, extracts its row from every `batch-inputs/*.txt` file, and runs:

```bash
python masquerade_interface.py \
  ${image} \            # 1: TIFF path
  ${spatial_metadata} \ # 2: spatial-annotation CSV
  ${outPath} \          # 3: output .ome.tiff path
  ${relevant_markers} \ # 4: marker metadata CSV, or "None"
  ${adjust_coords} \    # 5: True/False — crop image to coordinate bounds
  ${compress} \         # 6: True/False — bin/interpolate to target_size
  ${radius} \           # 7: circle radius in pixels
  ${filled} \           # 8: True/False — filled disk vs ring
  ${num_points} \       # 9: ring sample count (ignored when filled=True)
  ${preFilter_masks} \  # 10: True/False — spline filter before compression
  ${target_size}        # 11: target output size (GB)
```

See [masquerade_interface.py:13-34](masquerade/masquerade_interface.py) for the exact argument order.

After the Python step, the launcher locates the newly written `.ome.tiff` by its `basename`, `stat`s it, and — if ≥ 4 GB — runs `bfconvert` to produce the pyramidal variant.

#### batch-inputs/ format

Line-aligned, one row per sample across every file in [masquerade/batch-inputs/](masquerade/batch-inputs/):

| File | Contents |
|---|---|
| `qptiff-batch.txt` | absolute path to OME-TIFF / QPTIFF |
| `spatial_anno-batch.txt` | absolute path to spatial-annotation CSV (`x`, `y`, `cluster` columns) |
| `outPaths-batch.txt` | absolute path to target `.ome.tiff` output |
| `marker-metadata-batch.txt` | absolute path to marker metadata CSV, or the literal string `None` |

**To add a sample:** append one line to each file above.

---

### spatial-dynamics

#### What it does

Computes spatial relationships between annotated cell types. Two modules are selectable from the config:

- **`module=0` — Pairwise log-odds (PWLO).** For every ordered pair of cell types, counts how many cells of type B sit in an annulus with inner radius `p1` and outer radius `p2` microns around each cell of type A, then computes the log-odds of that co-occurrence against the background rate. Optionally adds a Kolmogorov–Smirnov effect-size test and a circos plot.
- **`module=1` — Multi-cell-type spatial circuit enrichment.** Takes a comma-delimited list of cell types (`circuit`, e.g. `pDC,CD4T,CD8T`) and computes n-way co-localization z-scores across neighborhoods, then optimizes neighborhood masks across a range of thresholds and emits an animation.

#### Inputs

- **Spatial-annotation CSV** — must contain `x`, `y`, and `cluster` columns (where `cluster` is the cell-type label).
- `resolution` config value — microns per pixel, used to convert `p1` / `p2` from microns to pixels.

#### Outputs

- **module=0:** per-sample CSVs with `(celltype_A, celltype_B, log_odds[, effect_size])` rows under each `${out_dir}`, plus optional circos plots when `draw=True`.
- **module=1:** neighborhood metadata, per-threshold masks, and a GIF generated by [optimize-neighborhoods.py](spatial-dynamics/optimize-neighborhoods.py).

#### Entry point and orchestration

`sbatch run-spatial_dynamics.s` reads [config-spatial_dynamics.txt](spatial-dynamics/config-spatial_dynamics.txt) and dispatches on `module`:

- `module=0` submits [run-pwlo.s](spatial-dynamics/run-pwlo.s) as an array job. Each task loads the `machinelearning` conda env, extracts its row, and runs:
  ```bash
  python3 run-pwlo.py \
    ${spatial_obj} \  # spatial-annotation CSV
    ${out_dir} \      # output directory
    ${label} \        # sample label
    ${resolution} \   # µm / pixel
    ${p1} ${p2}       # inner / outer annulus radius (µm)
  ```
  Note that `run-pwlo.py` currently hard-codes `draw=False` and `compute_effect_size=False` ([run-pwlo.py:13](spatial-dynamics/run-pwlo.py)) even though the config exposes those keys — edit the wrapper if you need them.

- `module=1` submits [run-spatial_circuit-enrichment.s](spatial-dynamics/run-spatial_circuit-enrichment.s), which calls:
  ```bash
  python3 run-spatial_circuit-enrichment.py \
    ${spatial_obj} ${out_dir} ${label} ${circuit}
  ```
  See **Known limitations** below — the circuit launcher does not extract per-array-task rows from the batch-inputs lists and the dispatch block in [run-spatial_dynamics.s:16-19](spatial-dynamics/run-spatial_dynamics.s) references a `${spatial_circuit_module_Path}` variable that is not defined in the shipped config.

#### batch-inputs/ format

Line-aligned, one row per sample across every file in [spatial-dynamics/batch-inputs/](spatial-dynamics/batch-inputs/):

| File | Contents |
|---|---|
| `spatial-annos.txt` | absolute path to spatial-annotation CSV |
| `outs.txt` | absolute path to output directory |
| `labels.txt` | sample label / ID |

**To add a sample:** append one line to each file above.

---

## Configuration reference

Every config file is just a shell-sourced `KEY=value` file — the launchers `source` it, so bash substitutions like `$(wc -l …)` work and absolute paths must be quoted if they contain spaces.

### [RunPhenomenalist/phenomenalist-config.txt](RunPhenomenalist/phenomenalist-config.txt)

```bash
configFile=phenomenalist-config.txt                 # self-reference, exported to child jobs

# Batch-input file paths (line-aligned; row N = sample N)
segmentation_file=batch-inputs/segmentation_files.txt  # per-sample segmentation CSV / TSV paths
failed_markers=batch-inputs/failed-markers.txt         # per-sample markers to drop (NULL if none)
nuclear_markers=batch-inputs/nuclear-markers.txt       # per-sample nuclear markers (NULL if none)
out_dir=batch-inputs/out_dirs.txt                      # per-sample output directories
classifier_labels=batch-inputs/labels.txt              # per-sample phenotyping labels (NULL if none)

# Phenotyping parameters
max_cells=3000000                                    # int — subsample threshold
HALO=T                                               # T/F — HALO vs Mesmer segmentation format
phenotyping_template=                                # optional gating template CSV (leave blank for none)
clustering_res=1,2                                   # comma list ('1,2,3') or range ('5:7' → c(5,6,7))

# SLURM parameters
batch_size=$(wc -l ${segmentation_file} | awk '{print $1}')   # auto-computed
module1_Path=run-phenomenalist.s                     # path to array-task launcher
module1_mem=300GB                                    # memory per array task
module1_time=0-12                                    # wall time (D-HH)
module1_partition=a100_short                         # SLURM partition (site-specific)
```

### [masquerade/configFile-batch.txt](masquerade/configFile-batch.txt)

```bash
configFile=configFile-batch.txt

# Batch-input file paths (absolute)
spatial_metadata=/path/to/batch-inputs/spatial_anno-batch.txt  # per-sample spatial-annotation CSV paths
image=/path/to/batch-inputs/qptiff-batch.txt                   # per-sample OME-TIFF / QPTIFF paths
outPath=/path/to/batch-inputs/outPaths-batch.txt               # per-sample output .ome.tiff paths
relevant_markers=/path/to/batch-inputs/marker-metadata-batch.txt  # per-sample marker metadata (or "None")

# Image preprocessing
adjust_coords=True        # True/False — crop image to coordinate bounds
compress=True             # True/False — bin/interpolate toward target_size
target_size=10            # int (GB) — target output size; triggers compression factor
preFilter_masks=True      # True/False — spline filter before compression

# Mask geometry
filled=True               # True/False — filled disk vs open ring
radius=5                  # int (pixels) — circle radius
num_points=50             # int — points sampled on ring (ignored when filled=True)

# SLURM parameters
sample_count=$(wc -l ${spatial_metadata} | awk '{print $1}')  # auto-computed
generate_masks_module_Path=run-masquerade-batch.sh            # path to array-task launcher
masquerade_partition=a100_short                               # SLURM partition (site-specific)
run_time=5                                                    # int (hours) — wall time
memory_=310GB                                                 # memory per array task
```

### [spatial-dynamics/config-spatial_dynamics.txt](spatial-dynamics/config-spatial_dynamics.txt)

```bash
configFile=config-spatial_dynamics.txt

# Batch-input file paths (line-aligned; row N = sample N)
spatial_annos=batch-inputs/spatial-annos.txt   # per-sample spatial-annotation CSV paths
outs=batch-inputs/outs.txt                     # per-sample output directories
labels=batch-inputs/labels.txt                 # per-sample sample IDs

# PWLO parameters
draw=False                 # True/False — emit circos plot (note: currently hard-coded False in run-pwlo.py)
resolution=0.3774          # float (µm/pixel) — instrument resolution
p1=3                       # float (µm) — inner annulus radius
p2=30                      # float (µm) — outer annulus radius
compute_effect_size=False  # True/False — run KS effect-size test (also hard-coded False in run-pwlo.py)

# Module selection
module=0                   # 0 = pairwise log-odds, 1 = multi-celltype neighborhoods
circuit=pdC,CD4T,CD8T      # comma-delimited target cell types (used when module=1)

# SLURM parameters
sample_count=$(wc -l ${spatial_annos} | awk '{print $1}')   # auto-computed
run_time=6                                                  # int (hours) — wall time
memory_=100GB                                               # memory per array task
spatial_dynamics_meta_partition=a100_short                  # SLURM partition (site-specific)
pwlo_module_Path=`pwd`/run-pwlo.s                           # path to PWLO array-task launcher
```

---

## Active vs. legacy files

The repository carries several older versions alongside the current implementations. Edit the **active** files; the legacy ones are retained for reference.

| Module | Active | Legacy / archived |
|---|---|---|
| RunPhenomenalist | `run-phenomenalist.R` + `RunPhenomenalist.R` | `v0/RunPhenomenalist-interface.R` (docopt-based) |
| masquerade (core) | `Masquerade.py` (class-based) | `Masquerade_v0.py` (function-based) |
| masquerade (CLI wrapper) | `masquerade_interface.py` | `masquerade_interface_v0.py`, `masquerade_interface_v1.py` (experimental — not wired into the launcher) |
| masquerade (helpers) | methods on the `Masquerade` class | `masquerade_utils.py` (superseded by `Masquerade.py`) |
| masquerade (launcher) | `run-masquerade-batch.sh` | `run-masquerade-batch_v0.sh` |
| spatial-dynamics | `pwlo_es_pt.py`, `n_simplex_neighborhoods.py`, `getNeighborhoods.py`, `optimize-neighborhoods.py` | `getNeighborhoods-v0.R`, `getNeighborhoods-v1.R`, `optimize-neighborhoods.R` |

---

## Run logs

Each module ships a `makeRunLog-batch.sh` that creates a timestamped file in `run-logs-batch/` capturing the config that was used plus the contents of every batch-input list. This is a post-run reproducibility record — not SLURM job metadata.

- [RunPhenomenalist/makeRunLog-batch.sh](RunPhenomenalist/makeRunLog-batch.sh) — invoke manually after a run; writes to [RunPhenomenalist/run-logs-batch/](RunPhenomenalist/run-logs-batch/).
- [masquerade/makeRunLog-batch.sh](masquerade/makeRunLog-batch.sh) — called automatically at the end of [run-masquerade-batch.sh](masquerade/run-masquerade-batch.sh:73); writes to [masquerade/run-logs-batch/](masquerade/run-logs-batch/).

SLURM's own `*_%j.err` / `*_%j.out` files land in the directory you ran `sbatch` from.

---

## Known limitations

- **Site-specific paths.** [configFile-batch.txt:3-9](masquerade/configFile-batch.txt) contains hard-coded `/gpfs/data/abl/tric/…` paths that must be edited before use elsewhere. The RunPhenomenalist CLI no longer needs editing for this — it discovers sibling files from its own script directory and reads the package source directory from `PHENOMENALIST_PKG_DIR` (with `library(phenomenalist)` and a historical GPFS path as fallbacks).
- **SLURM partitions are site-specific.** `a100_short` and `cpu_dev` will not exist on most clusters — edit each config and the `#SBATCH --partition=` line in [run-masquerade-batch.sh:3](masquerade/run-masquerade-batch.sh) before use.
- **Partial dependency manifests.** masquerade ships [environment.yml](masquerade/environment.yml) and RunPhenomenalist ships [environment.yml](RunPhenomenalist/environment.yml); spatial-dynamics does not. The [Prerequisites](#prerequisites) section lists every spatial-dynamics import so you can build its env yourself.
- **`run-pwlo.py` hard-codes `draw=False` and `compute_effect_size=False`** ([run-pwlo.py:13](spatial-dynamics/run-pwlo.py)) even though the config exposes those keys. Patch the wrapper if you need them.
- **`module=1` dispatch is incomplete.** [run-spatial_dynamics.s:17](spatial-dynamics/run-spatial_dynamics.s) submits `${spatial_circuit_module_Path}`, which is not defined in the shipped [config-spatial_dynamics.txt](spatial-dynamics/config-spatial_dynamics.txt). Additionally, [run-spatial_circuit-enrichment.s](spatial-dynamics/run-spatial_circuit-enrichment.s) does not extract `${spatial_obj}` / `${out_dir}` / `${label}` from the batch-input lists on a per-array-task basis. Expect to fix both before using the circuit-enrichment path.
- **`relevant_markers` is per-batch, not per-sample.** `run-masquerade-batch.sh` reads row N of `marker-metadata-batch.txt` like the other inputs, so to use one marker list across all samples you must repeat it on every line.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Array task N silently runs on the wrong sample | `batch-inputs/*.txt` files have mismatched line counts | `wc -l batch-inputs/*.txt` — every file should report the same number |
| masquerade fails at the `bfconvert` step | `java/17.0.0` not loaded, or the partition `#SBATCH` header overrides the env | Confirm `module load java/17.0.0` in [run-masquerade-batch.sh:8](masquerade/run-masquerade-batch.sh); run `java -version` inside an interactive SLURM session |
| `run-phenomenalist.R: error: RunPhenomenalist.R not found under …` | The wrapper could not auto-locate its siblings (rare — only happens when the script is copied without its directory, or run via a `source()` from another dir) | Set `PHENOMENALIST_DIR` to the directory containing `RunPhenomenalist.R` + `phenomenalist-utils.R` |
| `phenomenalist package not available` | Neither `PHENOMENALIST_PKG_DIR` is set nor `library(phenomenalist)` succeeds | Install the `phenomenalist` R package, or set `PHENOMENALIST_PKG_DIR` to its `R/` source directory |
| masquerade run succeeds but no `.ome.tiff` is found at the end | The launcher locates the output by running `ls -t ${out_dir}/${base}*.ome.tiff` — wrong `out_dir` permissions or an unexpected `basename` will make it miss | Check the directory of `outPath` in [configFile-batch.txt](masquerade/configFile-batch.txt) and confirm [masquerade_interface.py](masquerade/masquerade_interface.py) wrote the file; see [run-masquerade-batch.sh:34-39](masquerade/run-masquerade-batch.sh) |
| spatial-dynamics runs the wrong analysis | `module` set to `0` when you wanted `1`, or vice versa | Edit `module=` in [config-spatial_dynamics.txt](spatial-dynamics/config-spatial_dynamics.txt) |
| `sbatch` rejects the job immediately | Partition name is site-specific (`a100_short`, `cpu_dev`) | Update every `*_partition` key in the configs and every `#SBATCH --partition=` header to match your cluster |
