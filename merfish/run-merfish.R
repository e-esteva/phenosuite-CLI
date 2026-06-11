#!/usr/bin/env Rscript
#
# run-merfish.R
# =============
# CLI entry point for the MERFISH spatial-transcriptomics pipeline. Runs the
# full Import -> QC -> Normalize -> Reduce -> Cluster -> Spatial -> DE -> Export
# workflow on a single sample (one expression matrix + one metadata file).
#
# Invocation (named flags; the form run-merfish.s emits):
#
#   Rscript run-merfish.R \
#       --expression-file=PATH --metadata-file=PATH --out-dir=DIR \
#       --x-col=COL --y-col=COL \
#       [--sample-id=ID] [--transpose] \
#       [QC / processing / clustering / spatial / DE / export options]
#
# Run with --help for the complete option list. Designed to be sourced-free and
# standalone-runnable for single-sample testing outside SLURM.

## ---------------------------------------------------------------------------
## Usage / help
## ---------------------------------------------------------------------------

usage_text <- function() {
"Usage:
  Rscript run-merfish.R --expression-file=PATH --metadata-file=PATH \\
                        --out-dir=DIR --x-col=COL --y-col=COL [options]

Runs the headless MERFISH pipeline on a single sample. The expression matrix is
read as cells x genes (first column = cell ID); pass --transpose for genes x cells.
The metadata's first column is the cell ID used to align with the expression matrix.

Required:
  --expression-file=PATH   Cell x gene matrix (CSV / TSV / .gz). First col = cell ID.
  --metadata-file=PATH     Per-cell metadata (CSV / TSV / .gz). First col = cell ID.
  --out-dir=DIR            Output directory (created if missing).
  --x-col=COL              Metadata column holding the cell X coordinate.
  --y-col=COL              Metadata column holding the cell Y coordinate.

Sample / IO:
  --sample-id=ID           Output filename prefix.            [default: out-dir basename]
  --transpose              Treat the expression matrix as genes x cells.

QC thresholds (omit / set 'none' to disable a filter):
  --qc-min-counts=N        [default: 10]      --qc-max-counts=N        [default: 50000]
  --qc-min-genes=N         [default: 5]       --qc-max-genes=N         [default: 2000]
  --area-col=COL           Cell area/volume column (enables area + density filters)
  --qc-min-area=N          [default: none]    --qc-max-area=N          [default: none]
  --qc-min-density=N       [default: none]    --qc-max-density=N       [default: none]
  --negctrl-col=COL        Blank / negative-control count column
  --qc-max-negctrl-ratio=N [default: none]

Processing:
  --norm-method=M          lognorm|cp10k|cellvol|sct|none     [default: lognorm]
  --vol-col=COL            Cell volume column (for cellvol)
  --scale-method=M         zscore|center|none                 [default: zscore]
  --hvg-method=M           variance|vst|all                   [default: variance]
  --n-hvg=N                [default: 2000]
  --n-pcs=N                [default: 20]
  --umap-neighbors=N       [default: 15]      --umap-min-dist=F        [default: 0.3]

Clustering:
  --cluster-k=N            [default: 15]      --cluster-res=F          [default: 0.8]
  --leiden-objective=M     modularity|CPM                     [default: modularity]

Spatial / DE:
  --nhood-k=N              [default: 15]
  --svg-k=N                [default: 10]      --svg-n-top=N            [default: 200]

Export:
  --export-spe=BOOL        Write SpatialExperiment .rds        [default: TRUE]
  --export-figures=BOOL    Write PDF figure panels             [default: TRUE]
  --fig-width=F            [default: 8]  --fig-height=F [default: 6]  --fig-dpi=N [default: 300]
  --seed=N                 RNG seed                            [default: 42]
  -h, --help               Show this help and exit.

Sentinels that all mean 'not set' for a numeric option: NULL, null, NA, none, None, ''.

Environment variables:
  MERFISH_DIR              Directory holding RunMerfish.R and merfish-utils.R.
                           Defaults to the directory containing this script.

Prerequisites (provision the shipped conda env):
  conda env create -f environment.yml && conda activate runmerfish
"
}

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

die <- function(...) {
  message("run-merfish.R: error: ", ...)
  quit(save = "no", status = 1)
}

is_unset <- function(x) {
  if (is.null(x)) return(TRUE)
  if (length(x) == 0) return(TRUE)
  if (length(x) > 1) return(FALSE)
  if (is.na(x)) return(TRUE)
  trimws(as.character(x)) %in% c("", "NULL", "null", "NA", "none", "None", "NONE")
}

num_or_null <- function(x, flag) {
  if (is_unset(x)) return(NULL)
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) die(flag, " must be numeric (got '", x, "')")
  v
}

num_or_default <- function(x, default, flag) {
  v <- num_or_null(x, flag)
  if (is.null(v)) default else v
}

str_or_null <- function(x) if (is_unset(x)) NULL else as.character(x)

as_bool <- function(x, default = TRUE) {
  if (is_unset(x)) return(default)
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

## ---------------------------------------------------------------------------
## Script-directory discovery (for sourcing sibling files portably)
## ---------------------------------------------------------------------------

script_dir <- local({
  env_dir <- Sys.getenv("MERFISH_DIR", unset = "")
  if (nzchar(env_dir)) {
    if (!dir.exists(env_dir)) die("MERFISH_DIR is set but does not exist: ", env_dir)
    return(normalizePath(env_dir))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file_arg) == 1 && file.exists(file_arg)) return(dirname(normalizePath(file_arg)))
  getwd()
})

## ---------------------------------------------------------------------------
## Argument parsing
## ---------------------------------------------------------------------------

raw <- commandArgs(trailingOnly = TRUE)
if (length(raw) == 0 || any(raw %in% c("-h", "--help"))) {
  cat(usage_text())
  quit(save = "no", status = if (length(raw) == 0) 1 else 0)
}

# Defaults mirror the Shiny module's initial control values.
opts <- list(
  expression_file = NULL, metadata_file = NULL, out_dir = NULL,
  sample_id = NULL, transpose = "FALSE",
  x_col = NULL, y_col = NULL,
  qc_min_counts = "10", qc_max_counts = "50000",
  qc_min_genes = "5",  qc_max_genes = "2000",
  area_col = NULL, qc_min_area = NULL, qc_max_area = NULL,
  qc_min_density = NULL, qc_max_density = NULL,
  negctrl_col = NULL, qc_max_negctrl_ratio = NULL,
  norm_method = "lognorm", vol_col = NULL, scale_method = "zscore",
  hvg_method = "variance", n_hvg = "2000", n_pcs = "20",
  umap_neighbors = "15", umap_min_dist = "0.3",
  cluster_k = "15", cluster_res = "0.8", leiden_objective = "modularity",
  nhood_k = "15", svg_k = "10", svg_n_top = "200",
  export_spe = "TRUE", export_figures = "TRUE",
  fig_width = "8", fig_height = "6", fig_dpi = "300", seed = "42"
)

# Map of --flag -> opts slot. Boolean --transpose is a valueless switch.
flag_slot <- list(
  "--expression-file" = "expression_file", "--metadata-file" = "metadata_file",
  "--out-dir" = "out_dir", "--sample-id" = "sample_id", "--transpose" = "transpose",
  "--x-col" = "x_col", "--y-col" = "y_col",
  "--qc-min-counts" = "qc_min_counts", "--qc-max-counts" = "qc_max_counts",
  "--qc-min-genes" = "qc_min_genes", "--qc-max-genes" = "qc_max_genes",
  "--area-col" = "area_col", "--qc-min-area" = "qc_min_area", "--qc-max-area" = "qc_max_area",
  "--qc-min-density" = "qc_min_density", "--qc-max-density" = "qc_max_density",
  "--negctrl-col" = "negctrl_col", "--qc-max-negctrl-ratio" = "qc_max_negctrl_ratio",
  "--norm-method" = "norm_method", "--vol-col" = "vol_col", "--scale-method" = "scale_method",
  "--hvg-method" = "hvg_method", "--n-hvg" = "n_hvg", "--n-pcs" = "n_pcs",
  "--umap-neighbors" = "umap_neighbors", "--umap-min-dist" = "umap_min_dist",
  "--cluster-k" = "cluster_k", "--cluster-res" = "cluster_res",
  "--leiden-objective" = "leiden_objective",
  "--nhood-k" = "nhood_k", "--svg-k" = "svg_k", "--svg-n-top" = "svg_n_top",
  "--export-spe" = "export_spe", "--export-figures" = "export_figures",
  "--fig-width" = "fig_width", "--fig-height" = "fig_height", "--fig-dpi" = "fig_dpi",
  "--seed" = "seed"
)
bool_switches <- c("--transpose")

i <- 1L
while (i <= length(raw)) {
  tok <- raw[i]
  if (!grepl("^--", tok)) die("unexpected positional argument: '", tok, "' (see --help)")
  if (grepl("=", tok, fixed = TRUE)) {
    key <- sub("=.*$", "", tok); val <- sub("^[^=]*=", "", tok); i <- i + 1L
  } else if (tok %in% bool_switches) {
    key <- tok; val <- "TRUE"; i <- i + 1L
  } else {
    key <- tok
    if (i + 1L > length(raw)) die("missing value for ", key)
    val <- raw[i + 1L]; i <- i + 2L
  }
  slot <- flag_slot[[key]]
  if (is.null(slot)) die("unknown option: ", key, " (see --help)")
  opts[[slot]] <- val
}

## ---------------------------------------------------------------------------
## Normalize + validate
## ---------------------------------------------------------------------------

if (is_unset(opts$expression_file)) die("--expression-file is required")
if (is_unset(opts$metadata_file))   die("--metadata-file is required")
if (!file.exists(opts$expression_file)) die("expression file not found: ", opts$expression_file)
if (!file.exists(opts$metadata_file))   die("metadata file not found: ", opts$metadata_file)
if (is_unset(opts$x_col)) die("--x-col is required")
if (is_unset(opts$y_col)) die("--y-col is required")

if (is_unset(opts$out_dir)) {
  opts$out_dir <- getwd()
  message("run-merfish.R: no --out-dir given; using ", opts$out_dir)
}
dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)

sample_id <- if (is_unset(opts$sample_id)) basename(normalizePath(opts$out_dir)) else opts$sample_id

params <- list(
  expression_file = opts$expression_file,
  metadata_file   = opts$metadata_file,
  out_dir         = opts$out_dir,
  sample_id       = sample_id,
  transpose       = as_bool(opts$transpose, FALSE),
  x_col           = as.character(opts$x_col),
  y_col           = as.character(opts$y_col),

  qc_min_counts        = num_or_null(opts$qc_min_counts, "--qc-min-counts"),
  qc_max_counts        = num_or_null(opts$qc_max_counts, "--qc-max-counts"),
  qc_min_genes         = num_or_null(opts$qc_min_genes,  "--qc-min-genes"),
  qc_max_genes         = num_or_null(opts$qc_max_genes,  "--qc-max-genes"),
  area_col             = str_or_null(opts$area_col),
  qc_min_area          = num_or_null(opts$qc_min_area, "--qc-min-area"),
  qc_max_area          = num_or_null(opts$qc_max_area, "--qc-max-area"),
  qc_min_density       = num_or_null(opts$qc_min_density, "--qc-min-density"),
  qc_max_density       = num_or_null(opts$qc_max_density, "--qc-max-density"),
  negctrl_col          = str_or_null(opts$negctrl_col),
  qc_max_negctrl_ratio = num_or_null(opts$qc_max_negctrl_ratio, "--qc-max-negctrl-ratio"),

  norm_method     = as.character(opts$norm_method),
  vol_col         = str_or_null(opts$vol_col),
  scale_method    = as.character(opts$scale_method),
  hvg_method      = as.character(opts$hvg_method),
  n_hvg           = num_or_default(opts$n_hvg, 2000, "--n-hvg"),
  n_pcs           = num_or_default(opts$n_pcs, 20, "--n-pcs"),
  umap_neighbors  = num_or_default(opts$umap_neighbors, 15, "--umap-neighbors"),
  umap_min_dist   = num_or_default(opts$umap_min_dist, 0.3, "--umap-min-dist"),

  cluster_k        = num_or_default(opts$cluster_k, 15, "--cluster-k"),
  cluster_res      = num_or_default(opts$cluster_res, 0.8, "--cluster-res"),
  leiden_objective = as.character(opts$leiden_objective),

  nhood_k   = num_or_default(opts$nhood_k, 15, "--nhood-k"),
  svg_k     = num_or_default(opts$svg_k, 10, "--svg-k"),
  svg_n_top = num_or_default(opts$svg_n_top, 200, "--svg-n-top"),

  export_spe     = as_bool(opts$export_spe, TRUE),
  export_figures = as_bool(opts$export_figures, TRUE),
  fig_width      = num_or_default(opts$fig_width, 8, "--fig-width"),
  fig_height     = num_or_default(opts$fig_height, 6, "--fig-height"),
  fig_dpi        = num_or_default(opts$fig_dpi, 300, "--fig-dpi"),
  seed           = num_or_default(opts$seed, 42, "--seed")
)

if (!params$norm_method %in% c("lognorm", "cp10k", "cellvol", "sct", "none"))
  die("--norm-method must be one of lognorm|cp10k|cellvol|sct|none")
if (!params$scale_method %in% c("zscore", "center", "none"))
  die("--scale-method must be one of zscore|center|none")
if (!params$hvg_method %in% c("variance", "vst", "all"))
  die("--hvg-method must be one of variance|vst|all")
if (!params$leiden_objective %in% c("modularity", "CPM"))
  die("--leiden-objective must be modularity or CPM")

show <- function(label, value) {
  pretty <- if (is.null(value)) "<none>" else paste(value, collapse = ",")
  message(sprintf("  %-22s : %s", label, pretty))
}
message("== RunMerfish config ==")
show("sample_id",       params$sample_id)
show("expression_file", params$expression_file)
show("metadata_file",   params$metadata_file)
show("out_dir",         params$out_dir)
show("x_col / y_col",   paste(params$x_col, params$y_col, sep = " / "))
show("norm_method",     params$norm_method)
show("hvg_method",      paste0(params$hvg_method, " (n=", params$n_hvg, ")"))
show("n_pcs",           params$n_pcs)
show("cluster",         paste0("k=", params$cluster_k, " res=", params$cluster_res,
                               " ", params$leiden_objective))
show("export_spe",      params$export_spe)
message("=======================")

## ---------------------------------------------------------------------------
## Load deps + source pipeline (deferred so --help / arg errors stay clean)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages(suppressWarnings({
  library(Matrix); library(RANN); library(igraph); library(ggplot2)
}))

utils_file    <- file.path(script_dir, "merfish-utils.R")
pipeline_file <- file.path(script_dir, "RunMerfish.R")
if (!file.exists(utils_file))    die("merfish-utils.R not found under '", script_dir, "' (set MERFISH_DIR)")
if (!file.exists(pipeline_file)) die("RunMerfish.R not found under '", script_dir, "' (set MERFISH_DIR)")
Sys.setenv(MERFISH_DIR = script_dir)
source(utils_file)
source(pipeline_file)

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

RunMerfish(params)
