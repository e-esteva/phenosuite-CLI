# phenosuite-CLI

Command-line interface for the core PhenoSuite pipelines, designed to run on SLURM-based HPC clusters. The repository bundles eight independent spatial-biology modules that operate on multiplexed tissue imaging and spatial-transcriptomics data (CODEX, multiplexed IF, HALO / Mesmer / QuPath segmentations, MERFISH):

- **[segmentation](segmentation/)** — cell segmentation from multi-channel tissue images with four interchangeable routes: Segment Anything (SAM), DeepCell Mesmer, StarDist, and Cellpose (Python).
- **[RunPhenomenalist](RunPhenomenalist/)** — cellular scaling, dimensionality reduction, and clustering from single-cell segmentation tables (R).
- **[masquerade](masquerade/)** — circular cluster-mask generation for QuPath overlays on OME-TIFF / QPTIFF images (Python).
- **[spatial-dynamics](spatial-dynamics/)** — pairwise log-odds and multi-cell-type neighborhood enrichment for spatial cell–cell relationships (Python).
- **[merfish](merfish/)** — headless MERFISH spatial-transcriptomics pipeline: QC, normalization, clustering, neighborhood enrichment, spatially variable genes, and differential expression from cell × gene matrices (R).
- **[neighborhood_analysis](neighborhood_analysis/)** — KNN niche composition matrix, LOO stability sweep to select optimal neighbourhood count K₂, and MiniBatchKMeans assignment across multi-sample SpatialExperiment cohorts (R + Python).
- **[pcf](pcf/)** — inhomogeneous pair correlation functions between cell types from Vectra-format annotation tables: per-interaction PCF curves, normalized-PCF tables, and reference-cell-type interaction violins (R).
- **[gemma-phenotyper](gemma-phenotyper/)** — LLM cluster annotation with a Gemma model from HuggingFace weights, running offline on GPU nodes (R + Python).

The **segmentation** module's centroid-table output feeds directly into **RunPhenomenalist**, whose downstream `mask-inputs/` output in turn feeds **masquerade** — the three modules form one segmentation → phenotyping → visualization chain.

Each module is a self-contained `sbatch`-driven array job: a config file declares parameters, line-aligned text files in `batch-inputs/` declare the samples, and one `sbatch …` command fans the batch out across array tasks.

---

## Table of contents

1. [Repository layout](#repository-layout)
2. [Prerequisites](#prerequisites)
3. [Quickstart](#quickstart)
4. [Modules](#modules)
   - [segmentation](#segmentation)
   - [RunPhenomenalist](#runphenomenalist)
   - [masquerade](#masquerade)
   - [spatial-dynamics](#spatial-dynamics)
   - [merfish](#merfish)
   - [neighborhood_analysis](#neighborhood_analysis)
   - [pcf](#pcf)
   - [gemma-phenotyper](#gemma-phenotyper)
5. [Image preprocessing (QuPath crops for SAM)](#image-preprocessing-qupath-crops-for-sam)
6. [Configuration reference](#configuration-reference)
7. [Active vs. legacy files](#active-vs-legacy-files)
8. [Run logs](#run-logs)
9. [Known limitations](#known-limitations)
10. [Troubleshooting](#troubleshooting)
11. [License & Citation](#license--citation)

---

## Repository layout

```
phenosuite-CLI/
├── README.md
├── .gitignore
│
├── segmentation/                              # Python pipeline: SAM / Mesmer / StarDist / Cellpose routing
│   ├── run_segmentation.py                   #   CLI wrapper (named flags; standalone-runnable)
│   ├── segmentation_utils.py                 #   image I/O, channel resolution, per-route runners, exports
│   ├── segmentation-config.txt               #   batch config (edit this)
│   ├── segmentation-meta.s                   #   SLURM entry point — `sbatch` this
│   ├── run-segmentation.s                    #   SLURM array task (called by meta)
│   ├── makeRunLog-batch.sh                   #   captures config + input lists per run
│   ├── environment.yml                       #   conda env spec (python=3.10 + torch/tensorflow + 4 backends)
│   ├── batch-inputs/                         #   line-aligned per-image inputs
│   │   ├── images.txt
│   │   ├── out_dirs.txt
│   │   ├── labels.txt
│   │   ├── methods.txt
│   │   ├── nuclear_channels.txt
│   │   └── membrane_channels.txt
│   └── run-logs-batch/                       #   timestamped run logs (generated)
│
├── RunPhenomenalist/                         # R pipeline: scaling + dimensionality reduction + clustering
│   ├── RunPhenomenalist.R                    #   core pipeline function
│   ├── run-phenomenalist.R                   #   CLI wrapper (parses 9 positional args)
│   ├── export-anndata.R                      #   standalone spe.rds -> spe.h5ad converter
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
├── spatial-dynamics/                         # Python pipeline: PWLO + neighborhoods
│   ├── pwlo_es_pt.py                         #   pairwise log-odds implementation
│   ├── run-pwlo.py                           #   CLI wrapper for PWLO
│   ├── n_simplex_neighborhoods.py            #   n-way (3+ celltype) neighborhood enrichment
│   ├── getNeighborhoods.py                   #   wrapper orchestrating neighborhood analysis
│   ├── run-spatial_circuit-enrichment.py     #   CLI wrapper for circuit enrichment
│   ├── optimize-neighborhoods.py             #   threshold optimization + animation
│   ├── config-spatial_dynamics.txt           #   batch config (edit this)
│   ├── run-spatial_dynamics.s                #   SLURM entry point — `sbatch` this
│   ├── run-pwlo.s                            #   SLURM array task (module=0)
│   ├── run-spatial_circuit-enrichment.s      #   SLURM array task (module=1)
│   ├── batch-inputs/                         #   line-aligned per-sample inputs
│   │   ├── spatial-annos.txt
│   │   ├── outs.txt
│   │   └── labels.txt
│   └── *-v0.R, *-v1.R, optimize-neighborhoods.R   # legacy R implementations (archived)
│
├── merfish/                                  # R pipeline: MERFISH spatial transcriptomics
│   ├── RunMerfish.R                          #   core pipeline orchestrator
│   ├── run-merfish.R                         #   CLI wrapper (named flags; standalone-runnable)
│   ├── merfish-utils.R                       #   QC / normalize / PCA / UMAP / Leiden / spatial-stats / DE / figures
│   ├── merfish-config.txt                    #   batch config (edit this)
│   ├── merfish-meta.s                        #   SLURM entry point — `sbatch` this
│   ├── run-merfish.s                         #   SLURM array task (called by meta)
│   ├── environment.yml                       #   conda env spec (r-base=4.1 + deps)
│   ├── makeRunLog-batch.sh                   #   captures config + input lists per run
│   ├── batch-inputs/                         #   line-aligned per-sample inputs
│   │   ├── expression_files.txt
│   │   ├── metadata_files.txt
│   │   ├── out_dirs.txt
│   │   └── sample_ids.txt
│   └── run-logs-batch/                       #   timestamped run logs (generated)
│
├── neighborhood_analysis/                    # R pipeline: KNN niche + LOO sweep + neighbourhood assignment
│   ├── run-neighborhood_analysis.R           #   CLI wrapper (named flags + positional; standalone-runnable)
│   ├── neighborhood_analysis-utils.R         #   Python backend (sklearn KNN / MiniBatchKMeans) + R fallback + plot helpers
│   ├── neighborhood_analysis-config.txt      #   batch config (edit this)
│   ├── neighborhood_analysis-meta.s          #   SLURM entry point — `sbatch` this
│   ├── run-neighborhood_analysis.s           #   SLURM array task (called by meta)
│   ├── makeRunLog-batch.sh                   #   captures config + input lists per run
│   ├── environment.yml                       #   conda env spec (r-base=4.4 + sklearn/scipy/numpy)
│   ├── batch-inputs/                         #   line-aligned per-run inputs
│   │   ├── rds_files.txt                     #     .rds SpatialExperiment paths (; within a run = multi-sample)
│   │   ├── celltype_cols.txt
│   │   ├── out_dirs.txt
│   │   ├── labels.txt
│   │   └── condition_maps.txt
│   └── run-logs-batch/                       #   timestamped run logs (generated)
│
├── pcf/                                      # R pipeline: pair correlation functions between cell types
│   ├── run-pcf.R                             #   CLI wrapper (named flags; standalone-runnable)
│   ├── pcf-utils.R                           #   spatstat PCF estimators + curve tables + AUC violins
│   ├── pcf-config.txt                        #   batch config (edit this)
│   ├── pcf-meta.s                            #   SLURM entry point — `sbatch` this
│   ├── run-pcf.s                             #   SLURM array task (called by meta)
│   ├── makeRunLog-batch.sh                   #   captures config + input lists per run
│   ├── environment.yml                       #   conda env spec (r-base=4.4 + spatstat/ggpubr)
│   ├── batch-inputs/                         #   line-aligned per-run inputs
│   │   ├── vectra_files.txt                  #     dir of Vectra CSVs, ;-separated CSVs, or a .txt list
│   │   ├── out_dirs.txt
│   │   ├── labels.txt
│   │   ├── celltypes.txt
│   │   └── ref_celltypes.txt
│   └── run-logs-batch/                       #   timestamped run logs (generated)
│
└── gemma-phenotyper/                         # R pipeline: LLM cluster annotation (Gemma)
    ├── RunGemmaPhenotyper.R                  #   core pipeline function
    ├── run-gemma-phenotyper.R                #   CLI wrapper (named flags)
    ├── gemma-utils.R                         #   prompts, Ollama client, harmonisation
    ├── vectra-export.R                       #   Vectra-format CSV writer
    ├── loom-export.R                         #   .loom writer
    ├── infer.py                              #   transformers/PEFT backend (local checkpoints)
    ├── environment.yml                       #   conda env spec (R + Bioconductor)
    ├── gemma-phenotyper-config.txt           #   batch config (edit this)
    ├── gemma-phenotyper-meta.s               #   SLURM entry point — `sbatch` this
    ├── run-gemma-phenotyper.s                #   SLURM array task (called by meta)
    ├── makeRunLog-batch.sh                   #   captures config + input lists per run
    ├── batch-inputs/                         #   line-aligned per-sample inputs
    │   ├── spe_files.txt
    │   ├── out_dirs.txt
    │   ├── labels.txt
    │   ├── cluster_cols.txt
    │   └── tissues.txt
    └── run-logs-batch/                       #   timestamped run logs (generated)
```

---

## Prerequisites

Every module ships a conda `environment.yml` — RunPhenomenalist installs an R runtime, segmentation / masquerade / spatial-dynamics install Python runtimes. The sbatch launchers `source activate` each env by name, so create the env with the shipped file and the launchers work as-is.

### Cluster

- A **SLURM** cluster. Partition names in the shipped configs (`a100_short`, `cpu_dev`) are site-specific — edit them for your cluster before running.

### Conda storage location

`conda`/`mamba` put both environments and the package cache under `$HOME/.conda` (`envs/` and `pkgs/`) by default — easy to blow a small home-directory quota once the Bioconductor and torch stacks used here are involved. Redirect both to lab storage **once, before creating any env below**:

```bash
conda config --add envs_dirs /gpfs/data/abl/tric/<you>/conda/envs
conda config --add pkgs_dirs /gpfs/data/abl/tric/<you>/conda/pkgs
```

This writes to `~/.condarc` and changes *where* a named env is stored, not how you refer to it, so nothing else in this repo needs to change: every launcher activates envs by name (`source activate masquerade`, `conda activate runphenomenalist`, …), and name-based activation searches all registered `envs_dirs`. `conda env create -f environment.yml` and `conda create -n <name> …` (used for the separate `gemma-infer` Python env) both pick it up automatically.

Equivalent one-off, if you'd rather not touch `~/.condarc` (e.g. inside a single shell or job script):

```bash
export CONDA_ENVS_PATH=/gpfs/data/abl/tric/<you>/conda/envs
export CONDA_PKGS_DIRS=/gpfs/data/abl/tric/<you>/conda/pkgs
```

Already created an env under `$HOME` before setting this? Move it rather than reinstalling:

```bash
mv ~/.conda/envs/runphenomenalist /gpfs/data/abl/tric/<you>/conda/envs/
conda info --envs   # confirm it now resolves at the new path
```

`pip` (used inside `gemma-infer`) keeps its own cache under `$HOME/.cache/pip` regardless of the conda settings above — redirect it too if that matters for your quota:

```bash
export PIP_CACHE_DIR=/gpfs/data/abl/tric/<you>/.cache/pip
```

Two more caches live outside conda entirely and are covered where they're introduced rather than here: HuggingFace model weights for gemma-phenotyper (`HF_HOME`, see [Prerequisites > gemma-phenotyper](#gemma-phenotyper)) and the Python env `zellkonverter`/`basilisk` builds for RunPhenomenalist's AnnData export (see [Known limitations](#known-limitations)) — the latter installs under the R library path by default, not `$HOME/.conda`, so check `BASILISK_EXTERNAL_DIR` in your installed `basilisk` version's docs if it needs to move too.

### segmentation

Create the conda env from the shipped spec — the launcher expects it to be named `segmentation` and runs `source activate segmentation`:

```bash
cd segmentation/
conda env create -f environment.yml
```

[environment.yml](segmentation/environment.yml) installs `python=3.10`, `pytorch` + `torchvision`, and, via pip, all four backends: `cellpose`, `stardist` + `tensorflow`, `deepcell` + `tensorflow`, and `segment-anything`. That pulls in both PyTorch and TensorFlow, which is heavy and occasionally version-fussy on shared HPC modules — if you only need one or two routes, trim the `pip:` list to just those packages (see the comment block at the bottom of the file).

Also required at runtime:

- **A SAM checkpoint** for `--method=sam` — download one of the ViT checkpoints from the [Segment Anything repo](https://github.com/facebookresearch/segment-anything#model-checkpoints) and point `sam_checkpoint=` at it in the config (or pass `--sam-checkpoint=PATH` directly). Not required for the other three routes.
- **A GPU** is strongly recommended for all four routes but not required — `gpu=AUTO` in the config auto-detects CUDA via `torch.cuda.is_available()` and falls back to CPU silently.

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
- Optional AnnData export: `zellkonverter` (pulls in `basilisk`, which provisions its own Python env on first use — see [Known limitations](#known-limitations))

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

Create the conda env from the shipped spec — the launchers expect it to be named `spatial-dynamics` and run `source activate spatial-dynamics`:

```bash
cd spatial-dynamics/
conda env create -f environment.yml
# or, if the env already exists:
conda env update -f environment.yml --prune
```

The spec pins `python=3.10` and pulls `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `pillow` from `conda-forge`. See [spatial-dynamics/environment.yml](spatial-dynamics/environment.yml).

### merfish

Provision the shipped conda env (named `runmerfish`):

```bash
cd merfish/
conda env create -f environment.yml
conda activate runmerfish
```

What [environment.yml](merfish/environment.yml) pins:

- `r-base=4.1`
- Core analysis: `r-matrix`, `r-rann`, `r-igraph`, `r-rspectra` (spectral UMAP init; `cmdscale` fallback if absent), `r-data.table` (fast gz-aware reader; base `read.delim` fallback)
- Figures: `r-ggplot2`, `r-viridis`, `r-rcolorbrewer`, `r-pheatmap`
- Export / provenance: `r-jsonlite`, `bioconductor-spatialexperiment`, `bioconductor-singlecellexperiment`, `bioconductor-summarizedexperiment`, `bioconductor-s4vectors`

Unlike RunPhenomenalist, no external GitHub package is required — the pipeline ships entirely in `merfish/`. The `run-merfish.s` launcher runs `module load r/4.1.2` (site-specific); on a conda-only system, remove that line and `conda activate runmerfish` before `sbatch`.

### neighborhood_analysis

Provision the shipped conda env (named `neighborhoodr`):

```bash
cd neighborhood_analysis/
conda env create -f environment.yml
conda activate neighborhoodr
```

What [environment.yml](neighborhood_analysis/environment.yml) pins:

- `r-base=4.4` + `r-reticulate`
- CRAN: `r-dplyr`, `r-ggplot2`, `r-glue`, `r-jsonlite`, `r-rann`, `r-scales`, `r-tidyr`
- Bioconductor: `SpatialExperiment`, `SingleCellExperiment`, `SummarizedExperiment`, `MatrixGenerics`
- Python: `python=3.11`, `numpy<2`, `scipy`, `scikit-learn`, `psutil` (optional — enables memory logging)

The Python packages activate the high-performance backend automatically at runtime. If sklearn / scipy / numpy are absent, the pipeline falls back silently to RANN (R KNN) + base `kmeans`. The `run-neighborhood_analysis.s` launcher runs `module load r/4.4.0` (site-specific); on a conda-only system, remove that line and `conda activate neighborhoodr` before `sbatch`.
### pcf

Provision the shipped conda env (named `runpcf`):

```bash
cd pcf/
conda env create -f environment.yml
conda activate runpcf
```

What [environment.yml](pcf/environment.yml) pins:

- `r-base=4.4`
- Spatial statistics: `r-spatstat` (metapackage), `r-spatstat.geom`, `r-spatstat.explore` — the PCF estimators themselves
- Figures: `r-ggplot2`, `r-ggpubr` (`ggviolin()` + `stat_compare_means()`), `r-gridextra`
- Utilities: `r-glue`, `r-jsonlite`

No Python is required. The `pcf-v2` Shiny app reached spatstat the long way round — R → `reticulate` → `vectra_lib_v4.py` → `rpy2` → spatstat — so the maths always ran in R; this module calls spatstat directly and needs neither `reticulate` nor `rpy2`. The `run-pcf.s` launcher runs `module load r/4.4.0` (site-specific); on a conda-only system, remove that line and `conda activate runpcf` before `sbatch`.

### gemma-phenotyper

Create the conda env from the shipped spec — the launchers expect it to be named `gemma-phenotyper`:

```bash
cd gemma-phenotyper/
conda env create -f environment.yml
# or, if the env already exists:
conda env update -f environment.yml --prune
```

The spec installs an R runtime plus `r-curl`, `r-jsonlite`, the SpatialExperiment/SingleCellExperiment Bioconductor stack, and `r-hdf5r` (only needed for `--export-loom`). See [gemma-phenotyper/environment.yml](gemma-phenotyper/environment.yml).

**The `hf` backend — the default, and the only one that works on an offline compute node — needs a second Python env** with torch/transformers, pointed at by `gemma_python` in the config. It is deliberately kept out of `environment.yml` so the R and Python stacks can be provisioned independently and a CUDA build of torch can be chosen per cluster:

```bash
conda create -n gemma-infer python=3.10 && conda activate gemma-infer
pip install torch --index-url https://download.pytorch.org/whl/cu121 \
    "transformers>=4.50.0" "peft>=0.10.0" "accelerate>=0.28.0" "bitsandbytes>=0.43.0"
```

`transformers>=4.50.0` is the floor where Gemma 3 landed in a stable release. Swap the `--index-url` for `.../whl/cpu` if you have no GPU (expect it to be very slow for a 12B model).

**Pre-download the weights on a login node.** Gemma checkpoints on HuggingFace are gated (you must accept Google's license once, while logged in), and compute nodes are usually offline — the array tasks run with `HF_HUB_OFFLINE=1` and will not fetch anything:

```bash
huggingface-cli login
huggingface-cli download google/gemma-3-12b-it-qat-q4_0-unquantized \
    --local-dir /gpfs/data/myLab/models/gemma-3-12b-it-qat-q4_0-unquantized
```

Gemma 4 is Apache-2.0 and ungated, so it needs no `login` — but it does need `transformers>=5.10.0`:

```bash
hf download google/gemma-4-12B-it-qat-q4_0-unquantized \
    --local-dir /gpfs/data/myLab/models/gemma-4-12B-it-qat-q4_0-unquantized
```

**Set `HF_HOME` to project storage first.** It defaults to `$HOME/.cache/huggingface`, and most clusters quota `$HOME` well below a 12B checkpoint (`gemma-4-12B-it-qat-q4_0-unquantized` is a single unsharded 23.9 GB `model.safetensors`). `hf_xet` also stages chunks there during the transfer, so it matters even when you pass `--local-dir`:

```bash
export HF_HOME=/gpfs/data/myLab/.cache/huggingface
```

Set the same value as `hf_home=` in the config; the array tasks export it, so it applies on the compute nodes rather than only in the shell that ran `sbatch`.

Useful checks before committing to a 24 GB transfer:

```bash
hf download <repo> --local-dir <dest> --dry-run   # exact files + sizes, downloads nothing
hf env                                            # prints the resolved HF_HUB_CACHE
df -h <dest>                                      # confirm free space
hf cache verify <repo> --local-dir <dest>         # verify integrity afterwards
```

Point `model_dir` at whichever directory you downloaded. Browse available Gemma models — including newer generations as Google releases them — at [ai.google.dev/gemma](https://ai.google.dev/gemma).

**The optional `ollama` backend** needs no Python at all, only a reachable server — convenient off-cluster, rarely viable on one:

```bash
curl -fsSL https://ollama.com/install.sh | sh
ollama serve
ollama pull gemma3:12b-it-qat
```

---

## Quickstart

All modules follow the same pattern:

0. **One-time setup:** create the module's conda env — RunPhenomenalist and gemma-phenotyper also need a one-time package / model-weights step. If you're quota-constrained on `$HOME`, redirect conda's storage **before** running any `conda env create` below — see [Conda storage location](#conda-storage-location). Full picture for everything else (Java for masquerade, gated weight downloads for gemma-phenotyper, etc.) is in [Prerequisites](#prerequisites); the essentials are repeated below.
1. Edit the config file (`*-config*.txt`) in the module directory.
2. Add one line per sample to each file in `batch-inputs/` — **line N across every file = sample N**. The array size is auto-computed from the line count of the primary input list.
3. Submit the meta launcher with `sbatch`.

### segmentation

```bash
cd segmentation/
# 1. Edit segmentation-config.txt (default_method, model names, SLURM params)
vi segmentation-config.txt
# 2. Populate batch-inputs/ — one line per image in each file; use 'NULL' for empty rows
vi batch-inputs/images.txt
vi batch-inputs/out_dirs.txt
vi batch-inputs/labels.txt
vi batch-inputs/methods.txt
vi batch-inputs/nuclear_channels.txt
vi batch-inputs/membrane_channels.txt
# 3. Submit. segmentation-meta.s pre-validates that every batch-inputs file
#    has >= batch_size rows before calling sbatch.
sbatch segmentation-meta.s
```

To run the underlying CLI on a single image outside of SLURM (useful for debugging):

```bash
python run_segmentation.py --help
python run_segmentation.py \
    --image=/data/sample1.ome.tiff \
    --out-dir=/data/sample1/out-segmentation \
    --method=cellpose \
    --nuclear-channel=DAPI \
    --membrane-channel=CD45
```

### RunPhenomenalist

```bash
cd RunPhenomenalist/
# 0. One-time: create the conda env and install the phenomenalist package.
#    (The sbatch launcher itself runs `module load r/4.1.2`, not this env —
#    but the package must still be installed somewhere `library(phenomenalist)`
#    can find it, and this is the env the single-sample CLI below expects.)
conda env create -f environment.yml
conda activate runphenomenalist
Rscript -e 'remotes::install_github("igordot/phenomenalist")'
# 1. Edit phenomenalist-config.txt (clustering_res, SLURM params, etc.)
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
conda activate runphenomenalist
Rscript run-phenomenalist.R --help
Rscript run-phenomenalist.R \
    --segmentation-file=/data/sample1.csv \
    --out-dir=/data/sample1/out \
    --clustering-res=5:7 \
    --failed-markers=AF-Ak154,DAPI
```

### masquerade

```bash
cd masquerade/
# 0. One-time: create the conda env (Java 17 for bfconvert is loaded automatically
#    by the launcher via `module load java/17.0.0` — see Prerequisites)
conda env create -f environment.yml
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
# 0. One-time: create the conda env
conda env create -f environment.yml
# 1. Edit config-spatial_dynamics.txt — pick module=0 (PWLO) or module=1 (circuit)
vi config-spatial_dynamics.txt
# 2. Populate batch-inputs/ — one line per sample in each file
vi batch-inputs/spatial-annos.txt
vi batch-inputs/outs.txt
vi batch-inputs/labels.txt
# 3. Submit
sbatch run-spatial_dynamics.s
```

### merfish

```bash
cd merfish/
# 1. Edit merfish-config.txt (column mappings, QC/clustering params, SLURM params)
vi merfish-config.txt
# 2. Populate batch-inputs/ — one line per sample in each file; use 'NULL' for empty rows
vi batch-inputs/expression_files.txt
vi batch-inputs/metadata_files.txt
vi batch-inputs/out_dirs.txt
vi batch-inputs/sample_ids.txt
# 3. Submit. merfish-meta.s pre-validates that every batch-inputs file
#    has >= batch_size rows before calling sbatch.
sbatch merfish-meta.s
```

To run the underlying CLI on a single sample outside of SLURM (useful for debugging):

```bash
Rscript run-merfish.R --help
Rscript run-merfish.R \
    --expression-file=/data/sampleA/cell_by_gene.csv.gz \
    --metadata-file=/data/sampleA/cell_metadata.csv.gz \
    --out-dir=/data/results/sampleA \
    --sample-id=sampleA \
    --x-col=center_x --y-col=center_y \
    --area-col=cell_area --negctrl-col=blank_counts
```

### neighborhood_analysis

```bash
cd neighborhood_analysis/
# 1. Edit neighborhood_analysis-config.txt (K1, sweep range, SLURM params, etc.)
vi neighborhood_analysis-config.txt
# 2. Populate batch-inputs/ — one line per run; semicolons separate multiple RDS files within a run
vi batch-inputs/rds_files.txt       # e.g.: /data/s1.rds;/data/s2.rds;/data/s3.rds
vi batch-inputs/celltype_cols.txt   # column name, or NULL for auto-detect
vi batch-inputs/out_dirs.txt
vi batch-inputs/labels.txt
vi batch-inputs/condition_maps.txt  # e.g.: sample1=treated,sample2=control  (or NULL)
# 3. Submit. neighborhood_analysis-meta.s pre-validates row alignment before sbatch.
sbatch neighborhood_analysis-meta.s
```

To run on a single cohort outside of SLURM:

```bash
Rscript run-neighborhood_analysis.R --help

# count/pct mode — pooled random hold-out, no group structure assumed
Rscript run-neighborhood_analysis.R \
    --rds-files=/data/s1.rds;/data/s2.rds;/data/s3.rds \
    --out-dir=/data/results \
    --label=cohort1_20260620 \
    --celltype-col=celltype \
    --k1=10 \
    --loo-mode=count --loo-n=1 \
    --condition-map=s1=treated,s2=control,s3=control

# group mode — stratified hold-out (e.g. multiple timepoints, each with
# several replicates): holds out --loo-n samples from *every* condition
# group each fold, instead of a pooled random draw across all samples
Rscript run-neighborhood_analysis.R \
    --rds-files=/data/t0_a.rds;/data/t0_b.rds;/data/t0_c.rds;/data/t1_a.rds;/data/t1_b.rds;/data/t1_c.rds \
    --out-dir=/data/results \
    --label=timecourse_20260620 \
    --celltype-col=celltype \
    --k1=10 \
    --loo-mode=group --loo-n=1 \
    --condition-map=t0_a=T0,t0_b=T0,t0_c=T0,t1_a=T1,t1_b=T1,t1_c=T1
```

---

### pcf

```bash
cd pcf/
# 1. Edit pcf-config.txt (radius, instrument resolution, count threshold, SLURM params)
vi pcf-config.txt
# 2. Populate batch-inputs/ — one line per run; a run is one PCF analysis over one or more Vectra CSVs
vi batch-inputs/vectra_files.txt   # a directory of CSVs, /a.csv;/b.csv, or a .txt list
vi batch-inputs/out_dirs.txt
vi batch-inputs/labels.txt
vi batch-inputs/celltypes.txt      # comma-separated cell types, or NULL for the phenotypes shared by every file
vi batch-inputs/ref_celltypes.txt  # reference cell type for the interaction violins, or NULL
# 3. Submit. pcf-meta.s pre-validates row alignment before sbatch.
sbatch pcf-meta.s
```

To run on a single cohort outside of SLURM:

```bash
Rscript run-pcf.R --help

# every *.csv under a directory, all shared phenotypes, defaults elsewhere
Rscript run-pcf.R \
    --vectra-files=/data/vectra/cohort1 \
    --out-dir=/data/results/pcf \
    --label=cohort1_20260820 \
    --ref-celltype='CD8+ T cells'

# an explicit sample list and cell-type subset, 50 um radius on a 0.5 um/px scan
Rscript run-pcf.R \
    --vectra-files='/data/s1.csv;/data/s2.csv;/data/s3.csv' \
    --out-dir=/data/results/pcf \
    --label=cohort2_20260820 \
    --celltypes='CD8+ T cells,Macrophages,B cells' \
    --ref-celltype='CD8+ T cells' \
    --radius=50 --resolution=0.5 --count-threshold=10
```

---

### gemma-phenotyper

```bash
cd gemma-phenotyper/
# 0. One-time: create the R env, plus a separate Python env for the hf backend,
#    then download model weights on a login node (gated — see Prerequisites for
#    the full walkthrough, HF_HOME sizing, and the ollama alternative)
conda env create -f environment.yml
conda create -n gemma-infer python=3.10 && conda activate gemma-infer
pip install torch --index-url https://download.pytorch.org/whl/cu121 \
    "transformers>=4.50.0" "peft>=0.10.0" "accelerate>=0.28.0" "bitsandbytes>=0.43.0"
huggingface-cli login
huggingface-cli download google/gemma-3-12b-it-qat-q4_0-unquantized \
    --local-dir /gpfs/data/myLab/models/gemma-3-12b-it-qat-q4_0-unquantized
# 1. Edit gemma-phenotyper-config.txt — model_dir, prompt_style, GPU partition/gres
vi gemma-phenotyper-config.txt
# 2. Populate batch-inputs/ — one line per sample in each file; use 'NULL' for empty rows
vi batch-inputs/spe_files.txt
vi batch-inputs/out_dirs.txt
vi batch-inputs/labels.txt
vi batch-inputs/cluster_cols.txt
vi batch-inputs/tissues.txt
# 3. Submit. gemma-phenotyper-meta.s pre-validates row alignment AND the selected
#    backend (model tag set / checkpoint dir exists) before calling sbatch.
sbatch gemma-phenotyper-meta.s
```

Single sample, no SLURM:

```bash
Rscript run-gemma-phenotyper.R \
    --spe-file=/path/to/spe.rds \
    --out-dir=/path/to/out \
    --cluster-col=cluster_1 \
    --model-dir=/gpfs/data/myLab/models/gemma-3-12b-it-qat-q4_0-unquantized \
    --python=$(conda run -n gemma-infer which python) \
    --tissue='human lymph node'
```

---

## Modules

### segmentation

#### What it does

Segments cells from a multi-channel tissue image and routes the job through one of four interchangeable backends selected with `--method` (or the per-image `methods.txt` batch-input):

| Route | Backend | Typical use |
|---|---|---|
| `cellpose` (default) | [Cellpose](https://github.com/MouseLand/cellpose) | Generalist nuclei/cytoplasm segmentation; nuclear-only or two-channel "cyto" mode |
| `mesmer` | [DeepCell Mesmer](https://github.com/vanvalenlab/deepcell-tf) | Nuclear or whole-cell segmentation trained on multiplexed tissue images |
| `stardist` | [StarDist](https://github.com/stardist/stardist) | Star-convex nuclei segmentation; strong on dense, round nuclei |
| `sam` | [Segment Anything](https://github.com/facebookresearch/segment-anything) | Promptless automatic mask generation on a pseudo-RGB composite; no cell-type prior, boundaries are approximate |

The nuclear (and, for whole-cell routes, membrane) channel is resolved from the image by name or 0-based index, or auto-detected by pattern (`DAPI`, `Hoechst`, …) when left unset. All four routes emit the same output contract, so downstream tooling doesn't need to know which one produced a given mask.

**Memory scaling on large multiplex images.** QPTIFF / OME-TIFF panels (CODEX, multiplexed IF) routinely carry 20-60+ marker channels, each stored as its own TIFF page. `run_segmentation.py` reads channel *names* from metadata only, then decodes just the 1-2 pages it actually needs (nuclear ± membrane) — it never loads the full channel stack into memory. For a 40-channel QPTIFF, that's roughly a 20x reduction versus decoding every channel, which is what makes a 1-4 GB QPTIFF practical on a single array-task's memory budget (`module1_mem` in the config).

#### Inputs

- **Image** — multi-channel TIFF, OME-TIFF, or QPTIFF. Channel names are extracted from OME-XML metadata or the QPTIFF `Biomarker` tag, reusing the same detection logic as [masquerade](masquerade/Masquerade.py); plain TIFFs without embedded channel names fall back to `Channel_0`, `Channel_1`, … (pass `--nuclear-channel`/`--membrane-channel` as indices in that case).
- A SAM checkpoint `.pth` file, only for `--method=sam`.

#### Outputs

Written under `--out-dir`, prefixed with `--label` (defaults to the image basename):

| File | Contents |
|---|---|
| `{label}_mask.tif` | uint32 label mask, same H × W as the input image |
| `{label}_centroids.csv` | `label, y, x, area` — one row per segmented object. Having `label`/`y`/`x` columns present matches the bare **Mesmer** segmentation-table fingerprint that [RunPhenomenalist](#runphenomenalist) auto-detects, so this file (plus any marker intensities you join onto it) drops straight into RunPhenomenalist / masquerade regardless of which route produced it |
| `{label}_overlay.png` | nuclear channel with mask boundaries drawn in red, for quick QC |
| `{label}_provenance.json` | inputs, route, resolved channels, parameters, object count, timing, git SHA |

#### Entry point and orchestration

`sbatch segmentation-meta.s` reads [segmentation-config.txt](segmentation/segmentation-config.txt), validates that every `batch-inputs/*.txt` file has at least `batch_size` rows, and submits [run-segmentation.s](segmentation/run-segmentation.s) as an array job sized `1..batch_size` (with `--gres=${module1_gres}` attached when set). Each array task extracts its row from every `batch-inputs/*.txt` file with `sed -n "${SLURM_ARRAY_TASK_ID}p"`, resolves a `NULL`/empty `methods.txt` row to the config's `default_method`, and invokes the CLI with named flags:

```bash
python run_segmentation.py \
    --image="${image_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --label="${label_tmp}" \
    --method="${method_tmp}" \
    --nuclear-channel="${nuc_tmp}" \
    --membrane-channel="${mem_tmp}" \
    --diameter="${diameter}" \
    --cellpose-model="${cellpose_model}" \
    --stardist-model="${stardist_model}" \
    --sam-checkpoint="${sam_checkpoint}" \
    --sam-model-type="${sam_model_type}" \
    --resolution="${resolution}" \
    --min-size="${min_size}" \
    --tile-size="${tile_size}" \
    --tile-overlap="${tile_overlap}" \
    --seed="${seed}"
```

Full option reference:

```
Required:
  --image=PATH                 multi-channel TIFF / OME-TIFF / QPTIFF path
  --out-dir=DIR                output directory (created if missing)

Route:
  --method=STR                 sam | mesmer | stardist | cellpose         [default: cellpose]

Channels:
  --nuclear-channel=SPEC       channel name or 0-based index              [default: auto-detect]
  --membrane-channel=SPEC      channel name or 0-based index              [default: nuclear-only]

Route parameters:
  --diameter=FLOAT             cellpose expected object diameter, px      [default: 0 -> auto-estimate]
  --cellpose-model=STR         cellpose model_type                       [default: cyto3]
  --stardist-model=STR         pretrained StarDist2D model name           [default: 2D_versatile_fluo]
  --sam-checkpoint=PATH        SAM ViT checkpoint .pth (required for sam)
  --sam-model-type=STR         vit_h | vit_l | vit_b                     [default: vit_b]
  --sam-max-side=INT           refuse sam above N px on longest side     [default: 4096; 0 -> no limit]
  --resolution=FLOAT           microns/pixel, used by mesmer              [default: 0.5]

Post-processing:
  --min-size=INT               drop objects smaller than N px            [default: 0 -> no filter]
  --tile-size=INT               tile size, px (stardist/cellpose)         [default: 0 -> auto/whole-image]
  --tile-overlap=INT            tile overlap, px (stardist only)          [default: 64]

Misc:
  --label=STR                  output filename prefix                    [default: image basename]
  --gpu / --no-gpu              force GPU / CPU                          [default: auto-detect]
  --no-overlay                 skip the QC overlay PNG
  --no-centroids                skip the centroids CSV
  --seed=INT                    RNG seed                                 [default: 42]
  -h, --help                    show help and exit
```

Sentinels that mean *not set* for `--label`, `--nuclear-channel`, `--membrane-channel`, `--sam-checkpoint`, and per-image `batch-inputs/*.txt` rows: `NULL`, `null`, `NA`, `none`, `None`, `''` (empty). Unlike the other modules in this repo, **`0` is deliberately not a sentinel** here — it's a legitimate 0-based channel index (often the DAPI channel).

`run_segmentation.py` is fully standalone-runnable, so you can debug a single image on a login/dev node without SLURM. Route-specific backends (`torch`, `tensorflow`, `cellpose`, `stardist`, `deepcell`, `segment-anything`) are imported lazily inside each route's runner function, so `--help` and unrelated routes work even if only some backends are installed; a missing backend fails with an install hint rather than a raw traceback.

One environment variable pair controls provenance stamping:

| Variable | Meaning |
|---|---|
| `PHENOSUITE_GIT_SHA` | Recorded in `{label}_provenance.json`; `"unknown"` if unset |
| `PHENOSUITE_IMAGE_DIGEST` | Recorded in `{label}_provenance.json`; `"unknown"` if unset |

#### batch-inputs/ format

Every file in [segmentation/batch-inputs/](segmentation/batch-inputs/) holds **one row per image**, and row N must describe the same image across all files. Use the literal string `NULL` for a row with no value (never `0` — see above).

| File | Contents | Example row |
|---|---|---|
| `images.txt` | absolute path to the multi-channel image | `/data/sample1.ome.tiff` |
| `out_dirs.txt` | output directory (created if missing) | `/data/sample1/out-segmentation` |
| `labels.txt` | output filename prefix, or `NULL` (→ image basename) | `sample1` |
| `methods.txt` | `sam` \| `mesmer` \| `stardist` \| `cellpose`, or `NULL` (→ `default_method`) | `cellpose` |
| `nuclear_channels.txt` | channel name or index, or `NULL` (→ auto-detect) | `DAPI` |
| `membrane_channels.txt` | channel name or index, or `NULL` (→ nuclear-only) | `CD45` |

**To add an image:** append one line to each file above. The config's `batch_size=$(wc -l …)` will pick up the new count automatically.

---

### RunPhenomenalist

#### What it does

Reads single-cell segmentation measurements (one row per cell, one column per marker), builds a `SpatialExperiment` object, performs Leiden clustering at one or more resolutions, and emits marker-expression heatmaps, UMAP / spatial expression plots, and cluster assignments. Optionally assigns cell types using a manual gating template. The output `mask-inputs/` directory is formatted to feed directly into the **masquerade** module.

Supports HALO-, Mesmer-, and QuPath-style segmentation tables out of the box. The format is **auto-detected from the column headers** at runtime:

| Format | Fingerprint |
|---|---|
| **QuPath** | Space-delimited `Centroid X [µm\|px]` / `Centroid Y [µm\|px]` columns, and/or colon-separated per-compartment means (e.g. `CD3: Cell: Mean`). Centroid unit suffixes are tolerated — `µm`, `px`, `pixels`, or no suffix are all resolved to `x` / `y`. |
| **HALO** | `<marker> Cell/Nucleus/Cytoplasm/Membrane Intensity` columns and/or HALO metadata (`Classification`, `Completeness`). |
| **Mesmer** | Bare `label, y, x[, size, …]` header block followed by per-marker columns. |

Pass `--skip-cols=REGEX` to override the column-filter regex for custom schemas or when auto-detection picks wrong.

#### Inputs

- **Segmentation files** (CSV / TSV / TSV.gz): one row per cell, columns for `x`, `y`, and per-marker intensities. HALO, Mesmer, and QuPath layouts are all accepted — no flag needed.
- Optional **phenotyping template** CSV with manual gating rules (`phenotyping_template` key in the config).

#### Outputs

Written under each `${out_dir}/out-phenomenalist/`:

- `spe.rds` — serialized `SpatialExperiment` with expression, coordinates, and cluster assignments
- `spe.h5ad` — AnnData copy of the same object, written when `--export-anndata=true` (coordinates land in `adata.obsm['spatial']`); see [AnnData export](#anndata-export)
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
    --out-dir="${out_dir_tmp}" \
    --clustering-res="${clustering_res}" \
    --classifier-label="${classifier_label_tmp}" \
    --max-cells="${max_cells}" \
    --phenotyping-template="${phenotyping_template}" \
    --skip-cols="${skip_cols}" \
    --export-anndata="${export_anndata}"
```

Full option reference:

```
Required:
  --segmentation-file=PATH    per-cell segmentation CSV / TSV / TSV.gz
  --out-dir=DIR               output directory (created if missing)

Optional:
  --failed-markers=LIST       comma-separated markers to drop               [default: none]
  --nuclear-markers=LIST      comma-separated nuclear markers               [default: none]
  --classifier-label=LIST     comma-separated classifier labels             [default: none]
  --clustering-res=SPEC       '1,2,3' explicit list, or '5:7' range         [default: 1,2]
  --max-cells=N               cell subsample threshold                      [default: 100000]
  --phenotyping-template=PATH manual-gating template CSV                    [default: none]
  --skip-cols=REGEX           override column-filter regex (else auto)      [default: auto]
  --export-anndata=BOOL       also write spe.h5ad (AnnData) alongside spe.rds [default: false]
  -h, --help                  show help and exit
```

Sentinels that all mean *not set* for any option value: `NULL`, `null`, `NA`, `none`, `None`, `''` (empty), `0`. That makes every row of every `batch-inputs/*.txt` file safe to pad with `NULL`.

The legacy 8-positional invocation (`Rscript run-phenomenalist.R <seg> <failed> <nuclear> <out> <res> <label> <max> <template>`) still works for backward compatibility with any external caller, but the named-flag form is preferred.

The public [`phenomenalist`](https://github.com/igordot/phenomenalist) R package is loaded via `library(phenomenalist)` — no local sourcing of package `.R` files. If it's not installed the wrapper aborts with a clear install hint.

One environment variable controls local sibling-file discovery:

| Variable | Meaning |
|---|---|
| `PHENOMENALIST_DIR` | Directory holding `RunPhenomenalist.R` and `phenomenalist-utils.R`. Auto-detected from the script location; override only if you've split the files. |

`phenomenalist-meta.s` and `run-phenomenalist.s` both set `PHENOMENALIST_DIR` to their own directory automatically.

#### AnnData export

`spe.rds` can also be exported as an AnnData `.h5ad` file (via the Bioconductor [`zellkonverter`](https://bioconductor.org/packages/zellkonverter/) package) for use in scanpy / squidpy. `spatialCoords(spe)` is copied into a `reducedDim` named `"spatial"` before conversion, so it lands in `adata.obsm['spatial']`; the `counts` assay becomes `adata.X` and any other assays (`exprs`, `logcounts`, …) become `adata.layers`.

Two ways to get it:

1. **During the pipeline run** — pass `--export-anndata=true` (or set `export_anndata=True` in `phenomenalist-config.txt` for batch runs) and `spe.h5ad` is written next to `spe.rds` automatically.
2. **From an existing `spe.rds`**, without re-running clustering — use the standalone converter:
   ```bash
   Rscript export-anndata.R --spe-rds=/data/sample1/out/out-phenomenalist/sample1/spe.rds
   # writes spe.h5ad next to spe.rds; override with --out=PATH, and the exported
   # assay with --assay=NAME (default: counts)
   ```

Both paths call the same `export_anndata.mod()` helper in `phenomenalist-utils.R` and require the `zellkonverter` package — see [Known limitations](#known-limitations) for a note on first-use setup.

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

Reads a multi-channel tissue image (OME-TIFF or QPTIFF, routinely 4+ GB), draws a circular mask (filled disk or ring) around every cell centroid in a spatial-annotation CSV, and emits a multi-channel OME-TIFF containing the original biomarker channels as well as dedicated discrete channels holding masks for all clusters. The result drops straight into QuPath via its Bio-Formats reader for overlay visualization. Large outputs are automatically converted to tiled pyramidal OME-TIFF via the bundled `bfconvert`.

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

- `module=0` submits [run-pwlo.s](spatial-dynamics/run-pwlo.s) as an array job. Each task activates the `spatial-dynamics` conda env, extracts its row, and runs:
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

### merfish

#### What it does

Headless port of the PhenoSuite MERFISH module (`phenosuite/merfish/app.R`). Reads a per-sample cell × gene expression matrix and a matching cell-metadata table (with spatial coordinates) and runs the full segmentation-based spatial-transcriptomics workflow unattended, one sample per array task:

```
Import → QC → Normalize → HVG → PCA → UMAP → Leiden cluster →
Neighborhood enrichment → Spatially variable genes (Moran's I) → Cluster markers → Export
```

It emits cluster assignments, normalized expression, spatial statistics, figures, a full-provenance run manifest, and a `SpatialExperiment` `.rds` that drops straight into downstream PhenoSuite modules (phenotyping, PCF, spatial-interaction).

#### Inputs

- **Expression matrix** (CSV / TSV / `.gz`): cells × genes, first column = cell ID. Set `transpose=TRUE` in the config if your matrices are genes × cells.
- **Cell metadata** (CSV / TSV / `.gz`): first column = cell ID (used to align with the expression matrix), plus the X/Y coordinate columns named in the config and any optional area / negative-control / volume columns.

#### Outputs

Written under each `${out_dir}`, prefixed with `<sample_id>_`:

| File | Contents |
|---|---|
| `metadata_clusters.csv` | cell metadata + `cluster`, `Phenotype`, UMAP coordinates |
| `normalized_expression.csv` | normalized cell × gene matrix |
| `cluster_markers.csv` | one-vs-rest Wilcoxon markers per cluster |
| `svg.csv` | spatially variable genes (Moran's I, BH-adjusted p) |
| `neighborhood_enrichment.csv` | cluster × cluster co-localization Z-scores |
| `spe.rds` | `SpatialExperiment` (assays `exprs`/`logcounts`/`counts`, `UMAP`/`PCA` reducedDims, `Phenotype` colData) |
| `figures/*.pdf` | QC violin, spatial clusters, UMAP clusters |
| `run_manifest.json` | inputs, parameters, counts, timing, git SHA — full provenance |

#### Entry point and orchestration

`sbatch merfish-meta.s` reads [merfish-config.txt](merfish/merfish-config.txt), validates that every `batch-inputs/*.txt` file has at least `batch_size` rows (submission aborts with a clear error otherwise), and submits [run-merfish.s](merfish/run-merfish.s) as an array job sized `1..batch_size`. Each array task extracts its row from every list with `sed -n "${SLURM_ARRAY_TASK_ID}p"` and invokes the CLI with named flags built from the config:

```bash
Rscript run-merfish.R \
    --expression-file="${expression_file_tmp}" \
    --metadata-file="${metadata_file_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --sample-id="${sample_id_tmp}" \
    --x-col="${x_col}" --y-col="${y_col}" \
    --norm-method="${norm_method}" \
    --cluster-k="${cluster_k}" --cluster-res="${cluster_res}" \
    …  # full QC / processing / spatial / export options from the config
```

Key option groups (run `Rscript run-merfish.R --help` for the complete list):

```
Required:
  --expression-file=PATH   cell × gene matrix (first col = cell ID)
  --metadata-file=PATH     per-cell metadata (first col = cell ID)
  --out-dir=DIR            output directory (created if missing)
  --x-col=COL --y-col=COL  metadata coordinate columns

QC:        --qc-min/max-counts, --qc-min/max-genes, --area-col, --qc-min/max-area,
           --qc-min/max-density, --negctrl-col, --qc-max-negctrl-ratio
Process:   --norm-method (lognorm|cp10k|cellvol|sct|none), --vol-col, --scale-method,
           --hvg-method (variance|vst|all), --n-hvg, --n-pcs, --umap-neighbors, --umap-min-dist
Cluster:   --cluster-k, --cluster-res, --leiden-objective (modularity|CPM)
Spatial:   --nhood-k, --svg-k, --svg-n-top
Export:    --export-spe, --export-figures, --fig-width/height/dpi, --seed
Other:     --transpose (genes × cells input), --sample-id, -h/--help
```

Sentinels that all mean *not set* for any numeric option: `NULL`, `null`, `NA`, `none`, `None`, `''` (empty). That makes every row of every `batch-inputs/*.txt` file safe to pad with `NULL`. `run-merfish.R` is fully standalone-runnable, so you can debug a single sample on a login/dev node without SLURM.

One environment variable controls local sibling-file discovery:

| Variable | Meaning |
|---|---|
| `MERFISH_DIR` | Directory holding `RunMerfish.R` and `merfish-utils.R`. Auto-detected from the script location; override only if you've split the files. |

`merfish-meta.s` and `run-merfish.s` both set `MERFISH_DIR` to their own directory automatically.

#### batch-inputs/ format

Every file in [merfish/batch-inputs/](merfish/batch-inputs/) holds **one row per sample**, and row N must describe the same sample across all files. Use the literal string `NULL` for a row with no value.

| File | Contents | Example row |
|---|---|---|
| `expression_files.txt` | absolute path to cell × gene matrix | `/data/sampleA/cell_by_gene.csv.gz` |
| `metadata_files.txt` | absolute path to per-cell metadata | `/data/sampleA/cell_metadata.csv.gz` |
| `out_dirs.txt` | output directory (created if missing) | `/data/results/sampleA` |
| `sample_ids.txt` | output filename prefix, or `NULL` (→ out_dir basename) | `sampleA` |

The coordinate / area / negative-control / volume **column names** are global (set once in the config) and must exist in every sample's metadata file. `batch_size=$(wc -l < ${expression_file})` picks up the line count automatically.

---

### neighborhood_analysis

#### What it does

Headless port of the PhenoSuite NeighborhoodR module. Accepts a cohort of `SpatialExperiment` objects (multi-sample or single-sample), pools KNN niche composition vectors across all cells, runs an optional leave-one-out (LOO) stability sweep to select the optimal neighbourhood count K₂, then assigns every cell to a neighbourhood via MiniBatchKMeans and writes a single concatenated `SpatialExperiment` with per-cell `neighbourhood` and `sample` annotations.

```
Load SPEs → pool KNN niche matrix → LOO sweep (optional) →
MiniBatchKMeans assignment → cbind SPEs → export
```

Uses a high-performance scikit-learn / scipy backend (BallTree KNN, sparse scatter-add, MiniBatchKMeans). Falls back silently to RANN + base `kmeans` if sklearn / scipy / numpy are absent.

#### Inputs

- **SpatialExperiment `.rds` files** (one or more per run): each must contain spatial coordinates and a cell-type annotation column in `colData`. Multiple samples within one cohort are supplied as semicolon-separated paths on a single `batch-inputs/rds_files.txt` line.
- **Condition map** (optional): `sample1=treated,sample2=control` — maps sample names to experimental conditions for downstream comparative visualisation.

#### Outputs

Written under each `${out_dir}`, prefixed with `<label>_`:

| File | Contents |
|---|---|
| `<label>_joint_spe.rds` | Single concatenated `SpatialExperiment` with `neighbourhood` and `sample` in `colData` |
| `<label>_assignment_summary.csv` | Cell counts per sample × neighbourhood |
| `<label>_sweep_results.csv` | LOO stability scores for K₂ = `k2_min`…`k2_max` (empty when sweep is skipped) |
| `<label>_spatial_<sample>.png` | Per-sample tissue-space plot coloured by neighbourhood (when `make_plots=TRUE`) |
| `<label>_composition.png` | Stacked-bar neighbourhood composition per cell type (when `make_plots=TRUE`) |
| `<label>_provenance.json` | Inputs, parameters, R / Python versions, git SHA, image digest |

#### Entry point and orchestration

`sbatch neighborhood_analysis-meta.s` reads [neighborhood_analysis-config.txt](neighborhood_analysis/neighborhood_analysis-config.txt), validates that every `batch-inputs/*.txt` file has at least `batch_size` rows (submission aborts with a clear error otherwise), and submits [run-neighborhood_analysis.s](neighborhood_analysis/run-neighborhood_analysis.s) as an array job sized `1..batch_size`. Each array task extracts its row from every list and invokes the CLI:

```bash
Rscript run-neighborhood_analysis.R \
    --rds-files="${rds_files_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --label="${label_tmp}" \
    --celltype-col="${celltype_cols_tmp}" \
    --condition-map="${condition_maps_tmp}" \
    --k1="${k1}" --k2="${k2}" \
    --k2-min="${k2_min}" --k2-max="${k2_max}" \
    --loo-mode="${loo_mode}" --loo-n="${loo_n}" \
    --agg-fn="${agg_fn}" --condition-col="${condition_col}" \
    --seed="${seed}" $([ "${make_plots}" = "FALSE" ] && echo "--no-plots")
```

Full option reference:

```
Required:
  --rds-files=SPEC       semicolon-separated .rds paths, or path to a .txt file listing one path per line
  --out-dir=DIR          output directory (created if missing)
  --label=STR            output filename prefix

Optional — cell type:
  --celltype-col=SPEC    scalar column name, or per-sample map "s1=col1,s2=col2"  [default: auto-detect]

Optional — neighbourhood parameters:
  --k1=INT               K for KNN niche matrix                        [default: 10]
  --k2=INT               fixed neighbourhood count; skips LOO sweep    [default: NULL → run sweep]
  --k2-min=INT           LOO sweep lower bound                         [default: 3]
  --k2-max=INT           LOO sweep upper bound; NULL → k2_min+5        [default: NULL]

Optional — LOO sweep:
  --loo-mode=STR         'count', 'pct', or 'group' — see below                    [default: count]
  --loo-n=INT            hold-out count (meaning depends on --loo-mode)            [default: 1]
  --agg-fn=STR           'median' or 'mean' to aggregate fold scores               [default: median]

Optional — export:
  --condition-map=STR    "s1=cond1,s2=cond2" — condition labels for plots  [default: NULL]
  --condition-col=STR    colData column already holding condition labels    [default: NULL]
  --seed=INT             RNG seed                                           [default: 42]
  --no-plots             suppress PNG output
  -h, --help             show help and exit
```

`--loo-mode` controls how each fold's held-out sample set is drawn:

- **`count` / `pct`** — `--loo-n` samples (or percent of samples) drawn at random from the pooled sample list, ignoring any group structure.
- **`group`** — `--loo-n` samples drawn independently from *every* `--condition-col` / `--condition-map` group each fold (e.g. one replicate held out from every timepoint), rather than from the pooled list. Pooled random sampling doesn't reliably approximate this even when the total count is chosen to match: with N groups of R replicates, a random draw of size N lands as exactly one-per-group only a fraction of the time, and that fraction gets *worse* as N grows (≈16% for 4 groups of 3 replicates, ≈4% for 6 groups of 3) — the rest of the draws are uneven, sometimes wiping out a whole group while leaving another untouched. If the design has real replicate structure worth respecting, `group` mode is what actually tests it; `count`/`pct` only approximate it by luck. Requires `--condition-col` or every sample present in `--condition-map`.

Sentinels that all mean *not set*: `NULL`, `null`, `NA`, `none`, `None`, `''` (empty). Any row of any `batch-inputs/*.txt` file may safely use `NULL`. `run-neighborhood_analysis.R` is fully standalone-runnable without SLURM.

One environment variable controls sibling-file discovery:

| Variable | Meaning |
|---|---|
| `NEIGHBORHOODR_DIR` | Directory holding `neighborhood_analysis-utils.R`. Auto-detected from the script location; override only if the files are split across directories. |

#### batch-inputs/ format

Every file in [neighborhood_analysis/batch-inputs/](neighborhood_analysis/batch-inputs/) holds **one row per run** (where a "run" is one cohort of one or more samples). Row N must describe the same cohort across all files. Use `NULL` for any unused slot.

| File | Contents | Example row |
|---|---|---|
| `rds_files.txt` | semicolon-separated `.rds` SPE paths for this cohort | `/data/s1.rds;/data/s2.rds;/data/s3.rds` |
| `celltype_cols.txt` | column name in `colData`, `NULL` for auto-detect, or per-sample map | `celltype` |
| `out_dirs.txt` | output directory | `/data/results/cohort1` |
| `labels.txt` | output filename prefix | `cohort1_20260620` |
| `condition_maps.txt` | condition assignment string, or `NULL` | `s1=treated,s2=control,s3=control` |

**To add a cohort:** append one line to each file above.

---

### pcf

#### What it does

Headless port of the PhenoSuite `pcf-v2` module. Reads Vectra `cell_seg_data`-style annotation CSVs (one CSV = one sample), computes inhomogeneous pair correlation functions between every pair of cell types, and writes the PCF curves, a normalized-PCF table for one reference cell type, and the interaction plots the app produced.

```
Load Vectra CSVs → build marked point patterns → inhomogeneous PCF per
interaction (spatstat) → normalize by intensity → curve grid + AUC table +
interaction violins
```

For each sample the cells become a marked `ppp` (marks = phenotype) on the window spanned by the coordinates, and each pair of cell types is estimated with the matching spatstat estimator on a shared 500-point radius grid running from 0 to `radius / resolution` pixels:

| Pair | Estimator |
|---|---|
| `All` vs `All` | `pcfinhom()` — every cell against every cell |
| `All` vs celltype | `pcfdot.inhom()` — one type against the whole pattern |
| celltype vs celltype | `pcfcross.inhom()` — one type against another |

All three use `correction="isotropic"` with leave-one-out-free kernel intensities (`density(..., at="points", leaveoneout=FALSE)`), and each curve is divided by `Σ(1/λ₁)·Σ(1/λ₂) / (Δx·Δy)²` so curves are comparable across samples of different size and density — the normalization the app applied. A curve above 1 means the two types sit closer together than a random arrangement of the same intensities would put them.

The GUI ran R → `reticulate` → `vectra_lib_v4.py` → `rpy2` → spatstat, so every number it produced was already computed by spatstat in R. This module calls spatstat directly: same estimators, same arguments, same normalization, no Python in the loop.

#### Inputs

- **Vectra annotation CSVs** — comma- or tab-delimited, with `Cell X Position`, `Cell Y Position` and a phenotype column (`Phenotype` by default; `--phenotype-col` to change). `Tissue Category` is read when present and `Sample Name` is taken from the first row, falling back to the filename. **One CSV = one sample.** `gemma-phenotyper`'s `vectra_gemma_*.csv` export is written in exactly this format.
- **Cell types** (optional): `--celltypes` selects which phenotypes to analyse. The default is every phenotype present in *all* input files, matching the app's checkbox list.
- **Reference cell type** (optional): `--ref-celltype` is the type whose interactions are summarized in the AUC table and violins. Defaults to the first analysed cell type.

#### Outputs

Written under `${out_dir}/${label}/`:

| File | Contents |
|---|---|
| `<label>_pcf_summary.csv` | One row per sample × interaction: cell counts, `PCFsum`, normalization constant, and whether the pair was skipped |
| `<label>_pcf_curves.csv` | Long-format curves — one row per (sample, interaction, radius step) with raw `PCF` and `normPCF` |
| `<label>_ppc.rds` | Per-cell-type curve tables (the app's `ppc.rds`), one list element per cell type |
| `<label>-PCF_AUCs.csv` | `normPCF` series of every reference-cell-type interaction, one column per partner type plus `Sample` — the input format the PhenoSuite `pcf-builder` app reads for cross-cohort comparisons |
| `<label>-PCF-plots.pdf` | Curve grid: one panel per cell type, one line per partner, mean ± variance ribbon across samples |
| `<label>-PCF_AUC_violins.pdf` | Interaction violins vs the `All` baseline (single-sample runs) |
| `<label>-PCF_AUC_violins-global.pdf` | Same, pooled over every sample (multi-sample runs) |
| `<label>-PCF_AUC_violins-bySample.pdf` | Same, split by sample (multi-sample runs) |
| `individual-samples/<sample>-PCF_AUC_violins.pdf` | Per-sample violins (multi-sample runs) |
| `<label>_provenance.json` | Inputs, parameters, R / spatstat versions, git SHA, image digest |

#### Entry point and orchestration

`sbatch pcf-meta.s` reads [pcf-config.txt](pcf/pcf-config.txt), validates that every `batch-inputs/*.txt` file has at least `batch_size` rows and that `count_threshold >= 1` (submission aborts with a clear error otherwise), and submits [run-pcf.s](pcf/run-pcf.s) as an array job sized `1..batch_size`. Each array task extracts its row from every list and invokes the CLI:

```bash
Rscript run-pcf.R \
    --vectra-files="${vectra_files_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --label="${label_tmp}" \
    --celltypes="${celltypes_tmp}" \
    --ref-celltype="${ref_celltype_tmp}" \
    --radius="${radius}" --resolution="${resolution}" \
    --count-threshold="${count_threshold}" --min-count="${min_count}" \
    --phenotype-col="${phenotype_col}" \
    $([ "${make_plots}" = "FALSE" ] && echo "--no-plots")
```

Full option reference:

```
Required:
  --vectra-files=SPEC    directory of Vectra CSVs (searched recursively), ;-separated
                         CSV paths, or a .txt file listing one CSV path per line
  --out-dir=DIR          output directory (a --label subdirectory is created)

Optional — analysis:
  --label=STR            output filename prefix                    [default: YYYYMMDD]
  --celltypes=LIST       comma-separated cell types to analyse
                         [default: phenotypes present in every input file]
  --ref-celltype=NAME    reference cell type for the AUC violins   [default: first celltype]
  --radius=N             maximum radius in microns                 [default: 30]
  --resolution=N         instrument resolution, microns/pixel      [default: 0.377]
  --count-threshold=N    minimum cells of a type before its interactions
                         are computed (>= 1)                       [default: 10]
  --min-count=N          drop interactions whose smaller cell count is <= this
                         when assembling curves                    [default: 0]
  --phenotype-col=NAME   column holding the cell type labels       [default: Phenotype]

Optional — export:
  --no-plots             suppress PDF output
  -h, --help             show help and exit
```

`--radius` and `--resolution` set the radius grid together: curves run from 0 to `radius / resolution` pixels in 500 steps, and are plotted in microns. Getting `--resolution` wrong rescales every x-axis, so check it against the instrument that produced the scan (0.377 µm/px is the Vectra Polaris 20× default).

Interactions involving a cell type with fewer than `--count-threshold` cells in a given sample are skipped for that sample: the pair still appears in `<label>_pcf_summary.csv` with `skipped=TRUE` and no curve, and its column comes through as `NA` in that sample's `-PCF_AUCs.csv` rows. This is per sample, so a type that is abundant in one sample and rare in another contributes wherever it can.

Sentinels that all mean *not set*: `NULL`, `null`, `NA`, `none`, `None`, `''` (empty). Any row of any `batch-inputs/*.txt` file may safely use `NULL`. `run-pcf.R` is fully standalone-runnable without SLURM.

One environment variable controls sibling-file discovery:

| Variable | Meaning |
|---|---|
| `PCF_DIR` | Directory holding `pcf-utils.R`. Auto-detected from the script location; override only if the files are split across directories. |

#### Differences from the Shiny app

The estimators, arguments and normalization are unchanged. Four behaviours were fixed rather than reproduced, because reproducing them would have meant shipping known-wrong output:

- **AUC values are the raw `normPCF` series in every case.** The app's single-sample branch (its Python `pcf_AUC()`) multiplied them by `resolution` while its multi-sample branch did not, so one-sample and multi-sample runs came out in different units and could not be pooled.
- **AUC rows are assembled per sample**, so the `Sample` column always lines up with its values — including when a pair was skipped in one sample only. The app assumed every pair survived in every sample and would otherwise mislabel rows.
- **Curve-grid panels read their partner cell type per row** instead of inferring it from row position (same assumption), and each panel is titled with its own cell type instead of indexing the column names by panel number.
- **The dead `geom_signif()` layer was dropped.** Called without `comparisons=`, `stat_signif` has nothing to test and fails at draw time, so no bracket was ever drawn in the app either — only a warning per panel. The p-values that do appear come from `stat_compare_means(ref.group='All')`, which tests every interaction against the `All` baseline.

Also note that `-PCF_AUCs.csv` holds per-radius `normPCF` values, not integrated areas — "AUC" is the app's name for the file and is kept for `pcf-builder` compatibility.

#### pcf-builder compatibility

`<label>-PCF_AUCs.csv` is written to drop straight into the PhenoSuite `pcf-builder` app, which pools these tables across runs and plots one group against another. The format that app reads is:

| Requirement | How this module satisfies it |
|---|---|
| Loaded with `read.csv(…, row.names = 1)` | The table is written **with** row names, contiguous and unique |
| A `Sample` column | Multi-sample runs write the per-CSV sample name; single-sample runs write the run `--label`, so one run = one group |
| Every other column a cell type of numeric values | One column per partner cell type of the reference, plus `All` |
| Cell-type columns intersected across the loaded files | Column names are the phenotype labels themselves, so runs sharing a panel share columns |

Because that last step is an *intersection*, runs meant to be compared in `pcf-builder` should be produced with the same `--ref-celltype` and the same `--celltypes`. A reference cell type is excluded from its own AUC table, so two runs that used different references only share the columns neither one referenced — the app will silently offer the smaller list.

Two things are dropped before the file is written, because both are all-`NA` and would otherwise show up in the app as selectable-but-empty options: cell types whose interactions with the reference were skipped in *every* sample, and samples that contributed nothing at all. Each is reported on stderr when it happens, and both remain recorded in `<label>_pcf_summary.csv` with `skipped=TRUE`.

One case survives that filter by design: a cell type present in some samples and skipped in others keeps its real values and carries `NA` for the rest. In `pcf-builder`, selecting that cell type together with a sample that has no data for it leaves that group's violin empty and `stat_compare_means()` unable to use it as the reference group. Discarding those rows would mean discarding the samples where the interaction *was* measured, so they are kept.

Finally, `pcf-builder` renames samples into groups with `gsub(<sample name>, <new name>, Sample)` — it treats the `Sample` values as regular expressions and as substrings. Keep CSV filenames and `labels.txt` entries to letters, digits, `_` and `-`, and don't let one name be a prefix of another: `run+2` never matches itself, and renaming `s1` also rewrites part of `s10`. `run-pcf.R` prints a `NOTE` for both cases before writing the file, so a mis-grouped downstream plot is flagged at the point it becomes possible rather than discovered later.

#### batch-inputs/ format

Every file in [pcf/batch-inputs/](pcf/batch-inputs/) holds **one row per run** (where a "run" is one PCF analysis over one or more CSVs). Row N must describe the same run across all files. Use `NULL` for any unused slot.

| File | Contents | Example row |
|---|---|---|
| `vectra_files.txt` | directory of CSVs, `;`-separated CSV paths, or a `.txt` list | `/data/vectra/cohort1` |
| `out_dirs.txt` | output directory | `/data/results/pcf` |
| `labels.txt` | output filename prefix (and subdirectory name) | `cohort1_20260820` |
| `celltypes.txt` | comma-separated cell types, or `NULL` for the shared set | `CD8+ T cells,Macrophages,B cells` |
| `ref_celltypes.txt` | reference cell type, or `NULL` for the first one | `CD8+ T cells` |

**To add a run:** append one line to each file above.

---

### gemma-phenotyper

#### What it does

Assigns a cell-type label to every cluster (or every cell) in a `SpatialExperiment` / `SingleCellExperiment` by asking a Gemma model. It is the batch counterpart to PhenoSuite's `gemma_phenotyper` Shiny app, and shares its prompt-construction and inference code.

Two backends, selected with `backend=` in the config:

| Backend | How it runs | Needs |
|---|---|---|
| `hf` (default) | HuggingFace-format weights — a merged checkpoint or LoRA adapter, fine-tuned or an off-the-shelf base — loaded in-process by `infer.py`. Runs fully offline once the weights are on disk, so this is the HPC option. | `model_dir` + a Python env with torch/transformers (`gemma_python`), and a GPU for anything above ~4B |
| `ollama` | HTTP call per cluster to a running Ollama server. No weights loaded in the job. | `ollama_host` reachable **from the compute nodes** and a pulled `ollama_model` |

`local` is accepted as an alias for `hf`.

#### Prompt styles

Prompt style is independent of backend, set with `prompt_style=`:

| Style | What the model sees | Use when |
|---|---|---|
| `zero-shot` (default) | Tissue context + only the most distinctive markers for that cluster (top/bottom 10% by mean value, highest first) | The model was **not** fine-tuned on your ontology — e.g. an off-the-shelf `gemma-3-*-it-qat-q4_0-unquantized` base |
| `ontology` | Every marker as `marker=value`, plus an explicit ontology table | The checkpoint **was** fine-tuned against that ontology |

`zero-shot` sends the surviving markers highest-value-first, and states which assay the numbers came from (`exprs` arcsinh-scaled if present, else `scaled`, else z-scored `counts`) so the model isn't reading bare values with no sense of range:

> *System:* "You are an expert immunologist."
> *User:* "What cell type is described by `FOXP3:3.428, CD25:1.750, CD4:1.178, CD20:-0.147`? This is a `{tissue}`. Respond with ONLY JSON in the form {"cell_type": "\<3-word answer\>", "confidence": \<0-1\>} — no markdown, no explanation."

`ontology` sends every marker as `marker=value` alongside the ontology text:

> *System:* "You are a cellular phenotyping assistant for mIF data. Return JSON only."
> *User:* "Phenotype this cell cluster.\nCluster: `{id}` (n=`{n}` cells)\nMean markers (`{assay}`): `CD3=2.145, CD4=1.882, ...`\nOntology: `{ontology_text}`"

#### Model compatibility across Gemma generations

`infer.py` does not hardcode a Gemma generation. It reads the architecture class named in the checkpoint's own `config.architectures` and looks it up on the installed `transformers`, and decides multimodality structurally (does the config carry a `vision_config` / `audio_config` / `text_config`?) rather than by matching a `model_type` string. This matters because every generation so far has changed both:

| Checkpoint | `model_type` | Architecture class |
|---|---|---|
| `gemma-3-1b-it` | text-only | `Gemma3ForCausalLM` |
| `gemma-3-12b-it-qat-q4_0-unquantized` | `gemma3` | `Gemma3ForConditionalGeneration` |
| `gemma-4-E2B/E4B/31B` | `gemma4` | `Gemma4ForCausalLM` / `Gemma4ForConditionalGeneration` |
| `gemma-4-12B-it-qat-q4_0-unquantized` | `gemma4_unified` | `Gemma4UnifiedForConditionalGeneration` |

So a newer checkpoint generally needs **a new enough `transformers`, not a new pipeline**. If the installed `transformers` doesn't know the architecture, `infer.py` logs which class it wanted and falls back to an Auto class rather than crashing — but upgrading `transformers` is the real fix.

| Generation | `transformers` floor | Gating |
|---|---|---|
| Gemma 3 | `>=4.50.0` | Gated — accept Google's license, `huggingface-cli login` before downloading |
| Gemma 4 | `>=5.10.0` | Apache-2.0, ungated |

`gemma-4-12B`'s config declares `transformers_version: 5.10.0.dev0` — a **major** version jump from Gemma 3's floor, so treat that upgrade as a deliberate step (verify your other pinned deps) rather than a drop-in.

Multimodal checkpoints are loaded with their full architecture but prompted text-only — no images are ever sent.

**Gemma 4 thinking mode.** Gemma 4 models are reasoners with a configurable thinking mode that, when on, emits a `<|channel>thought …<channel|>` block *before* the answer — which would turn every prediction into a parse failure. `infer.py` pins `enable_thinking=False` in the chat template (falling back silently for Gemma 3 templates, which don't accept the kwarg), and additionally strips thinking blocks and ```` ``` ```` fences before parsing, then falls back to extracting the first balanced `{…}` span. So a stray reasoning trace or a "Sure, here's the JSON:" preamble still yields a usable label instead of an `unknown`.

#### Retrying failed clusters

A cluster is labelled `unknown` when it produced no usable answer — either a `request_failure` (the call errored or timed out) or a `parse_failure` (a reply came back but no `cell_type` could be extracted). Neither is a confidence judgement; there is no confidence threshold anywhere in the pipeline, and the `confidence: 0` on those rows is a hardcoded placeholder, not the model's own number.

Both are transient in principle, so the CLI retries them automatically: `--retry=N` (default 1) makes extra passes over just the failed prompts. On the `hf` backend the retry happens **inside the same `infer.py` process**, while the checkpoint is still resident — a fresh invocation would pay the full multi-minute load again to redo a handful of prompts. Set `--retry=0` to disable. Anything still failing afterwards is reported in the run log and keeps its `error`/`raw` columns in `*_gemma_predictions.csv`.

The Shiny app surfaces the same recovery as a **"Re-run Failed Clusters"** button in Step 4, which appears only when a run left failures. It re-requests just those prompts, merges any successes back into the predictions file, and then re-runs harmonisation, plots and exports so every downstream artefact reflects the recovered labels.

#### Marker glossary (strongly recommended)

A zero-shot model only knows markers it has read about. Anything named ambiguously, or cited mainly in specialist literature, it guesses at — **confidently**. Measured on a 39-plex human lymph-node panel with `gemma3:12b-it-qat`:

| Asked | Answered | Correct |
|---|---|---|
| which lineage is **TCF4**? | Regulatory T cells (0.95) | pDC |
| which lineage is **E2-2**? (same protein) | Plasma cells (0.95) | pDC |
| what co-expresses **TCF4 + IRF8 + CD123**? | Regulatory Tfh cells (0.95) | **pDC** |

TCF4 (E2-2) is the master transcription factor of plasmacytoid dendritic cells, and TCF4+IRF8+CD123 is a textbook pDC signature. The model reads "TCF4" as *T-cell factor* and commits at 0.95 — and renaming to E2-2 is differently wrong, not better. A pure 3,554-cell pDC cluster was consequently annotated *Macrophage*, recall **0.00**, in every prompt variant tested. No rewording fixes a wrong, confident prior.

Supplying the mapping does:

```bash
--marker-glossary=marker-glossary-example.txt     # CLI
marker_glossary=/path/to/glossary.txt             # config
```

On that same cluster: **`T follicular helper` (0.85) → `Plasmacytoid dendritic cell` (0.95)**.

**Write one for your own panel.** Copy [marker-glossary-example.txt](gemma-phenotyper/marker-glossary-example.txt) and edit. Format is `MARKER = short gloss`, one per line, `#` for comments; names match case- and punctuation-insensitively (`PD_1` / `PD-1` / `PD1` all hit). Only entries for a cluster's *elevated* markers are injected, so a long glossary costs no tokens on clusters that don't use those markers.

Include the markers a general-purpose model would plausibly get wrong:

- **transcription factors** — TCF4, IRF8, IRF4, PAX5, FOXP3
- **alias names that collide with other genes** — TCF4, S100, duplicate clones like CD103/CD103II
- **markers whose meaning is tissue-specific** — CD35, CD21, CD23 on follicular dendritic cells
- **anything you would argue about in a lab meeting**

Leave out markers with one unambiguous textbook meaning (CD3e, CD20, CD8) — the model already knows those and restating them only adds tokens.

#### Empirical confidence (N-sample voting)

The `confidence` the model returns is not informative — across a 49-cluster run it emitted ~0.95 for almost everything, including a cluster where every marker was below average. Set `vote_samples` in the config to replace it with a measured one:

```bash
vote_samples=10            # 0 disables
vote_temperature=0.7       # must be > 0; the main pass is greedy
vote_min_agreement=0.6     # fraction needed to replace the label
vote_override_floor=false
```

Each cluster that came back `Unclassified` is re-sampled `vote_samples` times, the modal label wins, and the agreement rate becomes the `confidence`. The vote distribution is written to `*_gemma_predictions.csv` as `vote_top`, `vote_agreement`, `vote_n`, `vote_basis` and `vote_detail`, so it is auditable even when the label doesn't change.

`vote_temperature` must exceed 0 — the main pass runs greedy (`temperature=0`) for reproducibility, so sampling at 0 would return N identical answers.

**Two kinds of `Unclassified` are treated differently.** When the `min_z` floor fired, the prompt literally says *no markers are elevated*; a majority label there is the model inventing a lineage from nothing, which is what the floor exists to prevent. Those are re-sampled for diagnosis (`vote_basis = no_elevated_markers`) but never overridden unless you set `vote_override_floor=true`. Clusters where the model declined *despite* having elevated markers (`vote_basis = model_declined`) can be overridden once agreement clears `vote_min_agreement`.

Voting runs **in-process on both backends** — inside `infer.py` for `hf`, so a 12B checkpoint is loaded once rather than once per pass. Cost scales directly: `vote_samples=10` means 10 extra generations per `Unclassified` cluster, so budget wall time accordingly on a large batch.

#### Harmonisation

Because each cluster is labelled in an independent call, the same population often comes back spelled several ways (`Regulatory T cell` / `Treg cell` / `T regulatory`). With `harmonize=true` (default) one extra generation maps every unique label to a canonical Title-Case form before results are written, and the mapping is saved to `*_harmonisation-map.csv`.

This runs on **both** backends. For `hf` it happens inside the same `infer.py` process, immediately after the main loop, so the checkpoint is loaded only once — a separate invocation would pay the full multi-minute load again just to map a handful of labels. If the model returns unparseable JSON the step degrades to keeping the raw labels rather than failing the run.

#### Inputs

| Input | Notes |
|---|---|
| `spe_files.txt` | `SpatialExperiment` or `SingleCellExperiment` `.rds`. Intensities are read from `exprs`, else `scaled`, else z-scored `counts`. |
| `cluster_cols.txt` | `colData` column holding cluster IDs. Required when `mode=cluster`. |
| `tissues.txt` | Tissue context injected into the prompt and used as the Vectra tissue fallback. Use `NULL` to omit. |

#### Outputs

Written to that sample's `out_dir`, prefixed with its label:

| File | Contents |
|---|---|
| `*_spe_gemma_annotated.rds` | input object plus `gemma_cell_type` / `gemma_confidence` in `colData` |
| `*_gemma_predictions.csv` | per-cluster model output, including `error`/`raw` for failed calls |
| `*_gemma_celltype-summary.csv` | cell counts per assigned type, descending |
| `*_gemma_colData.csv` | full `colData` dump |
| `*_harmonisation-map.csv` | raw → canonical label mapping (`ollama` + `harmonize=true`) |
| `vectra_gemma_*.csv` | Vectra-format export for PCF-toolkit compatibility (SpatialExperiment only) |
| `*_gemma_annotated.loom` | optional, with `export_loom=true` |
| `*_prompts.jsonl`, `*_predictions.jsonl` | exact prompts sent and raw responses — keep these for auditing |

#### Entry point and orchestration

[gemma-phenotyper-meta.s](gemma-phenotyper/gemma-phenotyper-meta.s) validates that every `batch-inputs/` file has at least `batch_size` rows **and** that the selected backend is usable, then submits [run-gemma-phenotyper.s](gemma-phenotyper/run-gemma-phenotyper.s) as `--array=1-${batch_size}`. Each array task extracts row `SLURM_ARRAY_TASK_ID` from every list and calls `run-gemma-phenotyper.R` with named flags. Only the flags the chosen backend reads are passed, so an unset value for the other backend can never be misread.

For `backend=hf` the pre-submit checks are: `model_dir` exists and holds a `config.json` or `adapter_config.json`; actual weight files (`*.safetensors` / `*.bin`) are present, since the tasks run offline and will not download them (a LoRA adapter is exempted, with a note that its base must already be cached); and `gemma_python` can import torch/transformers — a warning only, since compute nodes may load a different env. `--gres` is passed through from `module1_gres`; blank it for `ollama` so the job is not queued behind GPU availability it never uses.

For `backend=ollama`, reachability is checked from the *submit* host and is a warning, not a hard failure — compute nodes often have different network access, and that is what actually matters.

#### Running Ollama on a compute node (optional)

Prefer `backend=hf` on a cluster — it needs no network at all. If you nonetheless want the Ollama path and no shared server is reachable, start one inside the job and point the config at localhost. Add to the top of `run-gemma-phenotyper.s`:

```bash
ollama serve &
until curl -sf http://localhost:11434/api/tags >/dev/null; do sleep 2; done
ollama pull "${ollama_model}"
```

This costs a model load per array task, so it suits small batches; for larger ones prefer a persistent server on a node the cluster can reach.

#### batch-inputs/ format

Line-aligned, one row per sample across every file in [gemma-phenotyper/batch-inputs/](gemma-phenotyper/batch-inputs/):

| File | Contents |
|---|---|
| `spe_files.txt` | absolute path to the `.rds` object |
| `out_dirs.txt` | absolute path to output directory |
| `labels.txt` | sample label / output filename prefix |
| `cluster_cols.txt` | `colData` cluster column for that sample |
| `tissues.txt` | tissue context, or `NULL` |

**To add a sample:** append one line to each file above.

---

## Image preprocessing (QuPath crops for SAM)

`--method=sam` refuses images wider than `--sam-max-side` px on the longest side (default 4096 — see [Known limitations](#known-limitations)), because `SamAutomaticMaskGenerator` has no tiling and no cell-biology prior and doesn't scale to whole-slide QPTIFF/OME-TIFF images. `cellpose`, `mesmer`, and `stardist` don't need this — they tile or scale internally — so this section only matters if you specifically want to try the SAM route on a region of interest.

If your images are managed in QuPath (as with [masquerade](masquerade/)'s output), it can crop a region and export it as a channel-metadata-preserving OME-TIFF directly, which drops straight into `segmentation/run_segmentation.py --image=...`.

### Option A — GUI crop and export

1. Open the image in QuPath and draw a **rectangle annotation** around the region you want to segment (toolbar rectangle tool, or `R` shortcut).
2. Select the annotation, then **File → Export images...**
3. Choose **OME TIFF** as the export format — this preserves the per-channel `<Channel Name=...>` metadata that `segmentation`'s channel auto-detection and `--nuclear-channel`/`--membrane-channel` name lookup rely on. (A plain "TIFF" export typically drops channel names, forcing you to fall back to `--nuclear-channel` by index instead of name.)
4. Enable the region/crop option so the export is bounded to the selected annotation rather than the whole slide, and set downsample to `1` to keep full resolution.
5. Export, then point the segmentation module at the result:
   ```bash
   python run_segmentation.py --image=/path/to/cropped_roi.ome.tiff --out-dir=... --method=sam --sam-checkpoint=...
   ```

QuPath's exact dialog labels have shifted a little across versions (0.4.x vs 0.5.x) — if you don't see a region/crop toggle, check the current [QuPath documentation](https://qupath.readthedocs.io/) for your installed version.

### Option B — Groovy script (scripted / reproducible crops)

For a crop you want to reproduce exactly (e.g. as part of a batch of ROIs feeding `segmentation/batch-inputs/images.txt`), use QuPath's scripting console (**Automate → Show script editor**) or run headlessly via the `QuPath script` CLI. The building blocks are `RegionRequest` (defines the pixel bounding box) and `OMEPyramidWriter` (writes an OME-TIFF, optionally cropped to that region):

```groovy
import qupath.lib.regions.RegionRequest
import qupath.lib.images.writers.ome.OMEPyramidWriter

def server = getCurrentServer()
def x = 2000, y = 4000, w = 4000, h = 4000   // pixel bounding box of the ROI
def request = RegionRequest.createInstance(server.getPath(), 1.0, x, y, w, h)

new OMEPyramidWriter.Builder(server)
    .region(request)
    .tileSize(512)
    .build()
    .writeSeries('/path/to/cropped_roi.ome.tiff')
```

Treat this as a starting point rather than a drop-in for every QuPath version — the `OMEPyramidWriter` builder API has changed slightly across QuPath releases, so check the [scripting docs](https://qupath.readthedocs.io/) for your version if it doesn't run as-is. Keep `w`/`h` at or below `sam_max_side` (4096 by default) so the crop passes the segmentation module's guard without needing `--sam-max-side=0`.

---

## Configuration reference

Every config file is just a shell-sourced `KEY=value` file — the launchers `source` it, so bash substitutions like `$(wc -l …)` work and absolute paths must be quoted if they contain spaces.

### [segmentation/segmentation-config.txt](segmentation/segmentation-config.txt)

```bash
configFile=segmentation-config.txt

# Batch-input file paths (line-aligned; row N = image N)
image=batch-inputs/images.txt                       # multi-channel TIFF / OME-TIFF / QPTIFF path
out_dirs=batch-inputs/out_dirs.txt                   # output directory per image
labels=batch-inputs/labels.txt                       # output filename prefix (NULL -> image basename)
methods=batch-inputs/methods.txt                     # sam | mesmer | stardist | cellpose (NULL -> default_method)
nuclear_channels=batch-inputs/nuclear_channels.txt   # channel name or index (NULL -> auto-detect DAPI/Hoechst)
membrane_channels=batch-inputs/membrane_channels.txt # channel name or index (NULL -> nuclear-only)

# Route + shared parameters (global — apply to every image in the batch)
default_method=cellpose      # sam | mesmer | stardist | cellpose — used when methods.txt row is NULL
diameter=0                   # cellpose expected object diameter in px (0 = auto-estimate)
cellpose_model=cyto3
stardist_model=2D_versatile_fluo
sam_model_type=vit_b         # vit_h | vit_l | vit_b — must match sam_checkpoint
sam_checkpoint=              # path to a downloaded SAM ViT checkpoint .pth (required for method=sam)
sam_max_side=4096            # safety guard: refuse method=sam above this many px on the longest
                              # side (SAM has no tiling / cell-biology prior) — 0 disables the guard
resolution=0.5               # float (µm/pixel) — used by mesmer for cell sizing
min_size=15                  # drop objects smaller than this many px (0 = no filter)
tile_size=0                  # tile size in px for stardist/cellpose on large images (0 = auto/whole-image)
tile_overlap=64              # tile overlap in px (stardist only)
gpu=AUTO                     # TRUE | FALSE | AUTO — AUTO detects CUDA via torch
export_overlay=TRUE          # write a boundary-overlay PNG per image
export_centroids=TRUE        # write a label,y,x,area centroid CSV per image
seed=42

# SLURM parameters
batch_size=$(wc -l < ${image} | awk '{print $1}')
module1_Path=run-segmentation.s
module1_mem=64GB
module1_time=0-6
module1_partition=a100_short
module1_gres=gpu:1           # SLURM --gres value; leave blank for a CPU-only partition
```

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
phenotyping_template=                                # optional gating template CSV (leave blank for none)
skip_cols=                                           # optional regex override (leave blank for auto-detect)
clustering_res=1,2                                   # comma list ('1,2,3') or range ('5:7' → c(5,6,7))
export_anndata=False                                 # True/False — also write spe.h5ad (AnnData) alongside spe.rds

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

### [merfish/merfish-config.txt](merfish/merfish-config.txt)

```bash
configFile=merfish-config.txt

# Batch-input file paths (line-aligned; row N = sample N)
expression_file=batch-inputs/expression_files.txt   # per-sample cell × gene matrix paths
metadata_file=batch-inputs/metadata_files.txt       # per-sample cell-metadata paths
out_dir=batch-inputs/out_dirs.txt                    # per-sample output directories
sample_id=batch-inputs/sample_ids.txt               # per-sample output prefix (NULL → out_dir basename)

# Column mapping (global — must exist in every metadata file)
x_col=center_x                # cell X coordinate column
y_col=center_y                # cell Y coordinate column
area_col=cell_area            # cell area/volume column (NULL to disable area/density QC)
negctrl_col=blank_counts      # blank / negative-control count column (NULL to disable)
vol_col=cell_area             # cell volume column for cellvol normalization
transpose=FALSE               # TRUE if expression matrices are genes × cells

# QC thresholds (set any to NULL to disable that filter)
qc_min_counts=10              ; qc_max_counts=50000
qc_min_genes=5                ; qc_max_genes=2000
qc_min_area=NULL              ; qc_max_area=NULL
qc_min_density=NULL           ; qc_max_density=NULL
qc_max_negctrl_ratio=0.05

# Processing
norm_method=lognorm           # lognorm | cp10k | cellvol | sct | none
scale_method=zscore           # zscore | center | none
hvg_method=variance           # variance | vst | all
n_hvg=2000 ; n_pcs=20 ; umap_neighbors=15 ; umap_min_dist=0.3

# Clustering
cluster_k=15 ; cluster_res=0.8 ; leiden_objective=modularity   # modularity | CPM

# Spatial + DE
nhood_k=15 ; svg_k=10 ; svg_n_top=200

# Export
export_spe=TRUE ; export_figures=TRUE
fig_width=8 ; fig_height=6 ; fig_dpi=300 ; seed=42

# SLURM parameters
batch_size=$(wc -l < ${expression_file})   # auto-computed from the primary list
module1_Path=run-merfish.s                 # path to array-task launcher
module1_cpus=4                             # CPUs per array task
module1_mem=100GB                          # memory per array task
module1_time=0-06                          # wall time (D-HH)
module1_partition=a100_short               # SLURM partition (site-specific)
```

### [neighborhood_analysis/neighborhood_analysis-config.txt](neighborhood_analysis/neighborhood_analysis-config.txt)

```bash
configFile=neighborhood_analysis-config.txt

# Batch-input file paths (line-aligned; row N = cohort N)
rds_files=batch-inputs/rds_files.txt             # per-cohort RDS paths (; = multiple samples within a run)
celltype_cols=batch-inputs/celltype_cols.txt      # per-cohort celltype column name (NULL → auto-detect)
out_dirs=batch-inputs/out_dirs.txt               # per-cohort output directories
labels=batch-inputs/labels.txt                   # per-cohort output filename prefix
condition_maps=batch-inputs/condition_maps.txt   # per-cohort condition assignment (NULL if none)

# Neighbourhood parameters
k1=10           # int — K for KNN niche matrix
k2=NULL         # int — fixed neighbourhood count; NULL → run LOO sweep
k2_min=3        # int — LOO sweep lower bound (used when k2=NULL)
k2_max=NULL     # int — LOO sweep upper bound; NULL → k2_min+5

# LOO sweep options
loo_mode=count  # count | pct | group — how to partition samples in each fold
loo_n=1         # int — hold-out count, percentage, or per-group count (depends on loo_mode)
agg_fn=median   # median | mean — aggregate fold stability scores
condition_col=NULL  # colData column already holding condition labels (alternative to condition_map)

# Misc
seed=42
make_plots=TRUE

# SLURM parameters
batch_size=$(wc -l < ${rds_files} | awk '{print $1}')  # auto-computed from primary list
module1_Path=run-neighborhood_analysis.s               # path to array-task launcher
module1_mem=128GB                                       # memory per array task
module1_time=0-6                                        # wall time (D-HH)
module1_partition=cpu_short                             # SLURM partition (site-specific)
```

For a stratified sweep instead — e.g. holding out one replicate from every timepoint each fold, rather than a pooled random draw across all samples — change just these lines:

```bash
loo_mode=group
loo_n=1                  # samples held out per condition group, not a pooled total
condition_col=timepoint  # or leave NULL and give every sample a label in condition_maps.txt instead
```

`loo_mode=group` requires `condition_col` (or every sample present across `condition_maps.txt`) — `run-neighborhood_analysis.R` exits with a clear error at argument-parsing time if neither is set, and again mid-run if a sample ends up without a resolved label. `neighborhood_analysis-meta.s` doesn't check this itself — it only validates that batch-input files are row-aligned, not their contents.

---

### [pcf/pcf-config.txt](pcf/pcf-config.txt)

```bash
configFile=pcf-config.txt

# Batch-input file paths (line-aligned; row N = run N)
vectra_files=batch-inputs/vectra_files.txt    # per-run CSV directory, ;-separated CSVs, or .txt list
out_dirs=batch-inputs/out_dirs.txt            # per-run output directories
labels=batch-inputs/labels.txt                # per-run output prefix (also the subdirectory name)
celltypes=batch-inputs/celltypes.txt          # per-run cell types (NULL -> shared across all input files)
ref_celltypes=batch-inputs/ref_celltypes.txt  # per-run reference cell type (NULL -> first cell type)

# PCF parameters
radius=30                # float — maximum radius in microns
resolution=0.377         # float — instrument resolution, microns/pixel
count_threshold=10       # int   — minimum cells of a type before its interactions are computed (>= 1)
min_count=0              # int   — drop interactions whose smaller cell count is <= this
phenotype_col=Phenotype  # column holding the cell type labels
make_plots=TRUE          # write curve-grid + AUC violin PDFs

# SLURM parameters
batch_size=$(wc -l < ${vectra_files} | awk '{print $1}')  # auto-computed from primary list
module1_Path=run-pcf.s                                    # path to array-task launcher
module1_mem=64GB                                          # memory per array task
module1_time=0-6                                          # wall time (D-HH)
module1_partition=cpu_short                               # SLURM partition (site-specific)
```

Runtime scales with the number of interactions (`(n+1)(n+2)/2` per sample for `n` cell types) and with cell count per sample, since every pair fits its own kernel intensity estimate. Trim `celltypes.txt` to the populations you actually care about before raising `module1_time`.

---

### [gemma-phenotyper/gemma-phenotyper-config.txt](gemma-phenotyper/gemma-phenotyper-config.txt)

```bash
configFile=gemma-phenotyper-config.txt

# Batch-input file paths (line-aligned; row N = sample N)
spe_file=batch-inputs/spe_files.txt          # per-sample SpatialExperiment / SCE .rds paths
out_dir=batch-inputs/out_dirs.txt            # per-sample output directories
labels=batch-inputs/labels.txt               # per-sample output filename prefixes
cluster_cols=batch-inputs/cluster_cols.txt   # per-sample colData cluster column
tissues=batch-inputs/tissues.txt             # per-sample tissue context (or NULL)

# Backend: hf = HuggingFace weights in-process (offline-capable, the HPC option)
#          ollama = HTTP to a running server (needs compute-node reachability)
backend=hf

# Prompt style: zero-shot = tissue + most distinctive markers (not fine-tuned)
#               ontology  = all markers + ontology table (fine-tuned checkpoint)
prompt_style=zero-shot

# -- hf backend --
model_dir=/gpfs/data/myLab/models/gemma-3-12b-it-qat-q4_0-unquantized
base_model=google/gemma-3-1b-it      # only used when model_dir is a LoRA adapter
load_in_8bit=false                   # halves VRAM; requires bitsandbytes
ontology=                            # ontology text file (blank = built-in default)
gemma_python=python3                 # interpreter with torch/transformers/PEFT
hf_home=/gpfs/.../.cache/huggingface # HF cache root; exported into every array task

# -- ollama backend --
ollama_host=http://localhost:11434
ollama_model=gemma3:12b-it-qat
temperature=0.2

# Pipeline parameters (global)
mode=cluster               # cluster = one call per cluster; cell = one call per cell (slow)
markers=                   # comma-separated subset (blank = all rownames)
harmonize=true             # collapse near-duplicate labels
vote_samples=0             # N-sample vote on Unclassified clusters (0 = off, 10 typical)
vote_temperature=0.7       # must be > 0; main pass is greedy
vote_min_agreement=0.6     # fraction of votes needed to replace the label
vote_override_floor=false  # allow voting to override no-elevated-marker clusters
export_loom=false          # also write .loom alongside the annotated .rds

# SLURM parameters
batch_size=$(wc -l ${spe_file} | awk '{print $1}')   # auto-computed from primary list
module1_Path=run-gemma-phenotyper.s                  # array-task launcher
module1_mem=64GB                                     # host memory per task
module1_time=0-8                                     # wall time (D-HH)
module1_partition=gpu4_short                         # SLURM partition (site-specific)
module1_gres=gpu:1                                   # GPU request; blank for backend=ollama
```

Sizing for `backend=hf`: a 12B checkpoint in bf16 is roughly 24 GB, so the task needs a GPU with at least that much VRAM (A100 40/80GB, H100) or `load_in_8bit=true` to approximately halve it. Host memory should comfortably exceed the checkpoint size because the weights stream through RAM on load. For `backend=ollama` the job holds no weights at all — a CPU partition with ~16 GB and a blank `module1_gres` is plenty.

---

## Active vs. legacy files

The repository carries several older versions alongside the current implementations. Edit the **active** files; the legacy ones are retained for reference.

| Module | Active | Legacy / archived |
|---|---|---|
| segmentation | `run_segmentation.py` + `segmentation_utils.py` | — (no legacy files) |
| RunPhenomenalist | `run-phenomenalist.R` + `RunPhenomenalist.R` | `v0/RunPhenomenalist-interface.R` (docopt-based) |
| masquerade (core) | `Masquerade.py` (class-based) | `Masquerade_v0.py` (function-based) |
| masquerade (CLI wrapper) | `masquerade_interface.py` | `masquerade_interface_v0.py`, `masquerade_interface_v1.py` (experimental — not wired into the launcher) |
| masquerade (helpers) | methods on the `Masquerade` class | `masquerade_utils.py` (superseded by `Masquerade.py`) |
| masquerade (launcher) | `run-masquerade-batch.sh` | `run-masquerade-batch_v0.sh` |
| spatial-dynamics | `pwlo_es_pt.py`, `n_simplex_neighborhoods.py`, `getNeighborhoods.py`, `optimize-neighborhoods.py` | `getNeighborhoods-v0.R`, `getNeighborhoods-v1.R`, `optimize-neighborhoods.R` |
| neighborhood_analysis | `run-neighborhood_analysis.R` + `neighborhood_analysis-utils.R` | — (no legacy files) |
| pcf | `run-pcf.R` + `pcf-utils.R` | — (no legacy files) |

---

## Run logs

Each module ships a `makeRunLog-batch.sh` that creates a timestamped file in `run-logs-batch/` capturing the config that was used plus the contents of every batch-input list. This is a post-run reproducibility record — not SLURM job metadata.

- [segmentation/makeRunLog-batch.sh](segmentation/makeRunLog-batch.sh) — called automatically at the end of [run-segmentation.s](segmentation/run-segmentation.s); writes to `segmentation/run-logs-batch/`.
- [RunPhenomenalist/makeRunLog-batch.sh](RunPhenomenalist/makeRunLog-batch.sh) — invoke manually after a run; writes to [RunPhenomenalist/run-logs-batch/](RunPhenomenalist/run-logs-batch/).
- [masquerade/makeRunLog-batch.sh](masquerade/makeRunLog-batch.sh) — called automatically at the end of [run-masquerade-batch.sh](masquerade/run-masquerade-batch.sh:73); writes to [masquerade/run-logs-batch/](masquerade/run-logs-batch/).
- [merfish/makeRunLog-batch.sh](merfish/makeRunLog-batch.sh) — invoke manually before/after a run; writes to `merfish/run-logs-batch/`.
- [neighborhood_analysis/makeRunLog-batch.sh](neighborhood_analysis/makeRunLog-batch.sh) — called automatically at the end of [run-neighborhood_analysis.s](neighborhood_analysis/run-neighborhood_analysis.s); writes to `neighborhood_analysis/run-logs-batch/`.
- [pcf/makeRunLog-batch.sh](pcf/makeRunLog-batch.sh) — called automatically at the end of [run-pcf.s](pcf/run-pcf.s); writes to `pcf/run-logs-batch/`.
- [gemma-phenotyper/makeRunLog-batch.sh](gemma-phenotyper/makeRunLog-batch.sh) — invoke manually after a run; writes to [gemma-phenotyper/run-logs-batch/](gemma-phenotyper/run-logs-batch/). The per-sample `*_prompts.jsonl` / `*_predictions.jsonl` written into each `out_dir` are the finer-grained audit trail — they record the exact text sent to the model and its raw reply.

SLURM's own `*_%j.err` / `*_%j.out` files land in the directory you ran `sbatch` from.

---

## Known limitations

- **segmentation's `mesmer` and `sam` routes run whole-image inference only** — `--tile-size` / `--tile-overlap` are honored by `cellpose` and `stardist` (via `predict_instances_big`) but ignored (with a printed note) for `mesmer` and `sam`. Very large images may need to be cropped externally before running those two routes.
- **segmentation's `environment.yml` bundles all four backends into one env**, which pulls in both PyTorch and TensorFlow — heavy and occasionally version-fussy on shared HPC modules. Trim the `pip:` list to the route(s) you actually use if you hit install conflicts (see the comment block at the bottom of [segmentation/environment.yml](segmentation/environment.yml)).
- **segmentation's SAM route does not scale to whole-slide QPTIFFs, and is guarded accordingly.** `SamAutomaticMaskGenerator` has no tiling and no cell-biology prior — it runs Meta's ViT encoder once on a pseudo-RGB composite of the nuclear (+ membrane) channel. `run_sam()` refuses to run above `--sam-max-side` px on the longest side (default 4096) rather than silently burning GPU time on a slow, low-quality result; pass `--sam-max-side=0` (or a larger value) to override, but prefer `cellpose`, `mesmer`, or `stardist` for full QPTIFFs/OME-TIFFs and reserve `sam` for cropped ROIs.
- **Site-specific paths.** [configFile-batch.txt:3-9](masquerade/configFile-batch.txt) contains hard-coded `/gpfs/data/abl/tric/…` paths that must be edited before use elsewhere. The RunPhenomenalist CLI no longer needs editing for this — it discovers sibling files from its own script directory and loads `library(phenomenalist)` from the normal R library path.
- **SLURM partitions are site-specific.** `a100_short` and `cpu_dev` will not exist on most clusters — edit each config and the `#SBATCH --partition=` line in [run-masquerade-batch.sh:3](masquerade/run-masquerade-batch.sh) before use.
- **`run-pwlo.py` hard-codes `draw=False` and `compute_effect_size=False`** ([run-pwlo.py:13](spatial-dynamics/run-pwlo.py)) even though the config exposes those keys. Patch the wrapper if you need them.
- ** [run-spatial_circuit-enrichment.s](spatial-dynamics/run-spatial_circuit-enrichment.s) does not extract `${spatial_obj}` / `${out_dir}` / `${label}` from the batch-input lists on a per-array-task basis. Expect to fix both before using the circuit-enrichment path.
- **`relevant_markers` is per-batch, not per-sample.** `run-masquerade-batch.sh` reads row N of `marker-metadata-batch.txt` like the other inputs, so to use one marker list across all samples you must repeat it on every line.
- **AnnData export needs network access the first time it runs.** `zellkonverter` provisions its own isolated Python env via `basilisk` on first use. SLURM compute nodes are often offline, so do one `--export-anndata=true` run (or `Rscript export-anndata.R …`) on a login node with internet access before relying on the flag inside array jobs.
- **pcf curve panels are clipped to a PCF of 0–2.** The curve grid keeps the app's `ylim(0, 2)`, which *drops* points outside that band rather than zooming — a very strong short-range interaction leaves a visible gap at small radii instead of a spike. The unclipped values are in `<label>_pcf_curves.csv`.
- **pcf's `-PCF_AUCs.csv` holds per-radius `normPCF` values, not integrated areas.** The name is the Shiny app's and is kept so the files drop straight into `pcf-builder`; the violins therefore compare distributions of normalized PCF across the radius grid, not one AUC per sample.
- **gemma-phenotyper's `hf` backend reloads the checkpoint once per array task.** Each sample is an independent task, so a 12B checkpoint is loaded from disk for every one — several minutes of wall time per task before any inference starts. For a large batch of small samples that overhead can dominate; consider merging samples into fewer `.rds` objects, or raising `module1_time` accordingly.
- **`HF_HOME` defaults to `$HOME`, which is usually too small.** A 12B checkpoint is ~24 GB and cluster home directories are often quota'd below that. Set `hf_home=` in the config (it is exported into every array task) and export the same value in the shell you download from. A plain `model_dir` directory is read straight from disk and needs no cache — but a LoRA adapter's `base_model` is an HF repo id that must already be cached there, because the tasks run with `HF_HUB_OFFLINE=1`. `gemma-phenotyper-meta.s` warns pre-submit if that cache looks empty.
- **Gated weights must be pre-downloaded on a login node.** Gemma checkpoints require accepting Google's license while logged in, and the array tasks run with `HF_HUB_OFFLINE=1`. `gemma-phenotyper-meta.s` refuses to submit if `model_dir` has no `*.safetensors`/`*.bin`, but it cannot verify a LoRA adapter's *base* model is cached on the nodes — that one will only surface at runtime.
- **The `ollama` backend needs the server reachable from compute nodes, not just the submit host.** `gemma-phenotyper-meta.s` checks reachability from wherever you run `sbatch` and only warns on failure, because that check cannot speak for the compute nodes. If tasks come back with every cluster labelled `unknown` and `error=request_failure`, the nodes could not reach `ollama_host`.
- **A slow model can blow the per-request timeout (`ollama` backend).** Calls exceeding `--timeout` (default 600 s) are recorded as `request_failure` predictions labelled `unknown` rather than aborting the run — deliberate, so one slow cluster does not lose a whole sample, but a partially-`unknown` output is a real result you should check `*_gemma_predictions.csv` for.
- **`confidence` from the zero-shot backend is not calibrated.** In practice the model returns 0.95 for nearly every cluster regardless of evidence strength. Treat it as decorative, not as a quality filter.
- **Zero-shot accuracy is bounded by the marker panel.** Only the top/bottom 10% most distinctive markers per cluster reach the prompt, and populations that need markers outside that slice (or absent from the panel entirely) can be confidently mislabelled. Harmonisation makes labels *consistent*, not *correct* — spot-check against the marker values in `*_prompts.jsonl`.

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| Array task N silently runs on the wrong sample | `batch-inputs/*.txt` files have mismatched line counts | `wc -l batch-inputs/*.txt` — every file should report the same number |
| segmentation: `--method=sam requires 'torch' and 'segment-anything'` (or similar for cellpose/stardist/mesmer) | The route-specific backend isn't installed in the active env | Install just that backend (see the pip hint in the error), or `conda env create -f segmentation/environment.yml` for all four |
| segmentation: `--sam-checkpoint is required when --method=sam` | No checkpoint configured | Download a ViT checkpoint from the [Segment Anything repo](https://github.com/facebookresearch/segment-anything#model-checkpoints) and set `sam_checkpoint=` in the config (or pass `--sam-checkpoint=PATH`) |
| segmentation: `image is ... > --sam-max-side=4096` | SAM's whole-slide safety guard tripped — the image is larger than SAM can reasonably handle without tiling | Crop to a region of interest, switch to `cellpose`/`mesmer`/`stardist`, or override with `--sam-max-side=0` (no limit) / a larger value if you really want SAM on the full image |
| segmentation: `could not resolve a nuclear channel from [...]` | Auto-detection found no `DAPI`/`Hoechst`-like channel name (common on plain TIFFs with no embedded metadata) | Pass `--nuclear-channel` explicitly as a name or 0-based index — remember `0` is a valid index here, not a sentinel |
| segmentation produces very few / no objects | `--diameter` (cellpose) or the wrong channel resolved as nuclear | Check `{label}_overlay.png` and `{label}_provenance.json` for the resolved channel name; try `--diameter=0` for auto-estimation |
| masquerade fails at the `bfconvert` step | `java/17.0.0` not loaded, or the partition `#SBATCH` header overrides the env | Confirm `module load java/17.0.0` in [run-masquerade-batch.sh:8](masquerade/run-masquerade-batch.sh); run `java -version` inside an interactive SLURM session |
| `run-phenomenalist.R: error: RunPhenomenalist.R not found under …` | The wrapper could not auto-locate its siblings (rare — only happens when the script is copied without its directory, or run via a `source()` from another dir) | Set `PHENOMENALIST_DIR` to the directory containing `RunPhenomenalist.R` + `phenomenalist-utils.R` |
| `phenomenalist package not available` | Neither `PHENOMENALIST_PKG_DIR` is set nor `library(phenomenalist)` succeeds | Install the `phenomenalist` R package, or set `PHENOMENALIST_PKG_DIR` to its `R/` source directory |
| masquerade run succeeds but no `.ome.tiff` is found at the end | The launcher locates the output by running `ls -t ${out_dir}/${base}*.ome.tiff` — wrong `out_dir` permissions or an unexpected `basename` will make it miss | Check the directory of `outPath` in [configFile-batch.txt](masquerade/configFile-batch.txt) and confirm [masquerade_interface.py](masquerade/masquerade_interface.py) wrote the file; see [run-masquerade-batch.sh:34-39](masquerade/run-masquerade-batch.sh) |
| pcf: `no phenotype is present in all N input files` | The phenotype labels differ across CSVs (e.g. one sample was annotated with a different cluster naming) | Pass `--celltypes` explicitly, or harmonise the `Phenotype` column across exports |
| pcf: `no PCF curves were computed for reference cell type '…'` | Every interaction the reference type takes part in fell below `--count-threshold` in every sample | Lower `count_threshold`, or pick a more abundant `ref_celltype` |
| pcf-builder: a group's violin is empty, or `Can't find specified reference group` | The selected cell type was skipped (below `count_threshold`) in that sample, so its rows are `NA` | Pick another reference group, or re-run that cohort with a lower `count_threshold`; `<label>_pcf_summary.csv` lists which pairs were skipped |
| pcf-builder: renaming a sample into a group does nothing, or renames the wrong one | `Sample` values contain regex metacharacters, or one name is a prefix of another — the app renames with `gsub()` | Rename the input CSVs (or `labels.txt` entries) to letters/digits/`_`/`-`, with no name a prefix of another; `run-pcf.R` prints a `NOTE` when it writes such names |
| pcf: curves look shifted along the x-axis | `resolution` does not match the instrument that produced the scan | Set `resolution` (microns/pixel) to the acquisition value — it converts `radius` into the pixel grid the coordinates live on |
| spatial-dynamics runs the wrong analysis | `module` set to `0` when you wanted `1`, or vice versa | Edit `module=` in [config-spatial_dynamics.txt](spatial-dynamics/config-spatial_dynamics.txt) |
| `sbatch` rejects the job immediately | Partition name is site-specific (`a100_short`, `cpu_dev`) | Update every `*_partition` key in the configs and every `#SBATCH --partition=` header to match your cluster |
| `--export-anndata=true` / `export-anndata.R` errors with `'zellkonverter' package is required` | Package not installed in the active R library | `Rscript -e 'BiocManager::install("zellkonverter")'`, or re-provision the conda env from the updated [environment.yml](RunPhenomenalist/environment.yml) |
| AnnData export hangs or fails on a download/network error inside a SLURM job | `basilisk` (zellkonverter's Python backend) is trying to build its env on an offline compute node | Run once with `--export-anndata=true` on a login node with internet access to build the env first, then resubmit |

---

## License & Citation

If you use phenosuite-CLI in your research, please cite it:

> Esteva, E. (2026). phenosuite-CLI: Command-line interface for PhenoSuite spatial omics pipelines [Computer software]. GitHub. https://github.com/e-esteva/phenosuite-CLI

A machine-readable citation is available in [`CITATION.cff`](CITATION.cff). GitHub will automatically render a "Cite this repository" button from it.
