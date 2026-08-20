#!/usr/bin/env Rscript
#
# run-neighborhood_analysis.R
# ===========================
# CLI entry point for the NeighborhoodR spatial neighborhood analysis pipeline.
#
# Supports two invocation styles:
#
#   (1) Named flags (preferred):
#       Rscript run-neighborhood_analysis.R \
#           --rds-files=path1.rds;path2.rds \
#           --out-dir=DIR \
#           [--label=LABEL] \
#           [--celltype-col=COLNAME] \
#           [--k1=10] \
#           [--k2=N | --k2-min=3 --k2-max=N] \
#           [--loo-mode=count] [--loo-n=1] \
#           [--agg-fn=median] \
#           [--condition-col=COLNAME] \
#           [--condition-map=s1=cond1,s2=cond2] \
#           [--seed=42] \
#           [--no-plots]
#
#   (2) Positional (batch script back-compat):
#       Rscript run-neighborhood_analysis.R \
#           <rds_files> <celltype_col> <out_dir> <label> \
#           <k1> <k2_min> <k2_max> <loo_mode> <loo_n> \
#           <condition_map> <seed>
#
# --rds-files accepts:
#   • Semicolon-separated paths:  path1.rds;path2.rds;path3.rds
#   • Path to a .txt file with one path per line (when value ends in .txt)
#
# --condition-map accepts comma-separated sample=condition pairs:
#   sample1=treated,sample2=control
#
# Run with --help for full details.

## ---------------------------------------------------------------------------
## Usage / help
## ---------------------------------------------------------------------------

usage_text <- function() {
"Usage:
  Rscript run-neighborhood_analysis.R --rds-files=PATHS --out-dir=DIR [options]
  Rscript run-neighborhood_analysis.R <rds_files> <ct_col> <out_dir> <label>
      <k1> <k2_min> <k2_max> <loo_mode> <loo_n> <cond_map> <seed>

Builds a KNN niche composition matrix, runs a LOO stability sweep to select
the optimal neighborhood count K2, assigns cells with MiniBatchKMeans, and
saves a concatenated SpatialExperiment RDS with 'neighborhood' in colData.

Required:
  --rds-files=LIST      Semicolon-separated .rds SpatialExperiment paths, OR
                        path to a .txt file listing one path per line.
  --out-dir=DIR         Root output directory (a label subdirectory is created).

Options:
  --label=STR           Output file prefix.                  [default: YYYYMMDD]
  --celltype-col=NAME   Shared celltype column name, OR per-sample pairs:
                        's1=celltype,s2=cluster'. Omit for auto-detection.
  --k1=N                KNN neighbours for niche matrix.     [default: 10]
  --k2=N                Use this fixed K2 — skips LOO sweep. [default: auto]
  --k2-min=N            LOO sweep lower bound.               [default: 3]
  --k2-max=N            LOO sweep upper bound.               [default: n_celltypes-1]
  --loo-mode=STR        'count' or 'pct'.                    [default: count]
  --loo-n=N             Hold-out count or percentage.        [default: 1]
  --agg-fn=STR          'median' or 'mean' stability metric. [default: median]
  --condition-col=NAME  colData column with condition labels. [default: none]
  --condition-map=STR   Manual labels: 'sample1=ctrl,sample2=treated'.
  --seed=N              Random seed.                         [default: 42]
  --no-plots            Skip PNG plot generation.
  -h, --help            Show this help and exit.

Sentinels meaning 'not set': NULL, null, NA, none, None, '' (empty).

Environment variables:
  NEIGHBORHOODR_DIR    Directory of this script + utils (auto-detected).
  RETICULATE_PYTHON    Python binary (e.g. /opt/venv/bin/python).

Output (written to --out-dir/--label/):
  {label}_joint_spe.rds            Concatenated SPE with neighborhood in colData
  {label}_assignment_summary.csv   Sample x neighborhood cell counts
  {label}_sweep_results.csv        LOO stability scores per K2 (if sweep ran)
  {label}_spatial_{sample}.png     Per-sample spatial projection (unless --no-plots)
  {label}_composition_barplot.png  Celltype x neighborhood barplot (unless --no-plots)
  {label}_provenance.json          Full run provenance

Prerequisites:
  conda env create -f environment.yml && conda activate neighborhoodr
"
}

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

die <- function(...) {
  message("run-neighborhood_analysis.R: error: ", ...)
  quit(save="no", status=1)
}

is_unset <- function(x) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  if (length(x) > 1)                return(FALSE)
  if (is.na(x))                     return(TRUE)
  trimws(as.character(x)) %in%
    c("", "0", "NULL", "null", "NA", "na", "none", "None", "NONE")
}

parse_list <- function(x, sep=",") {
  if (is_unset(x)) return(NULL)
  parts <- trimws(unlist(strsplit(as.character(x), sep, fixed=TRUE)))
  parts[nzchar(parts)]
}

parse_int <- function(x, flag) {
  v <- suppressWarnings(as.integer(x))
  if (is.na(v)) die(flag, " must be an integer (got '", x, "')")
  v
}

parse_num <- function(x, flag) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) die(flag, " must be numeric (got '", x, "')")
  v
}

parse_rds_files <- function(x) {
  if (is_unset(x)) return(NULL)
  x <- trimws(as.character(x))
  if (grepl("\\.txt$", x, ignore.case=TRUE) && file.exists(x)) {
    paths <- trimws(readLines(x, warn=FALSE))
    return(paths[nzchar(paths)])
  }
  parts <- trimws(unlist(strsplit(x, ";", fixed=TRUE)))
  parts[nzchar(parts)]
}

parse_condition_map <- function(x) {
  if (is_unset(x)) return(NULL)
  pairs <- parse_list(x, sep=",")
  if (is.null(pairs)) return(NULL)
  result <- list()
  for (p in pairs) {
    kv <- strsplit(p, "=", fixed=TRUE)[[1]]
    if (length(kv) != 2)
      die("--condition-map: bad pair '", p, "' (expected sample=condition)")
    result[[trimws(kv[1])]] <- trimws(kv[2])
  }
  result
}

parse_ct_cols <- function(x) {
  if (is_unset(x)) return(NULL)
  if (!grepl("=", x, fixed=TRUE)) return(trimws(x))
  pairs  <- parse_list(x, sep=",")
  result <- character(0)
  for (p in pairs) {
    kv <- strsplit(p, "=", fixed=TRUE)[[1]]
    if (length(kv) != 2)
      die("--celltype-col: bad pair '", p, "' (expected sample=col)")
    result[trimws(kv[1])] <- trimws(kv[2])
  }
  result
}

## ---------------------------------------------------------------------------
## Script-directory discovery
## ---------------------------------------------------------------------------

script_dir <- local({
  env_dir <- Sys.getenv("NEIGHBORHOODR_DIR", unset="")
  if (nzchar(env_dir)) {
    if (!dir.exists(env_dir)) die("NEIGHBORHOODR_DIR does not exist: ", env_dir)
    return(normalizePath(env_dir))
  }
  args     <- commandArgs(trailingOnly=FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value=TRUE))
  if (length(file_arg) == 1 && file.exists(file_arg))
    return(dirname(normalizePath(file_arg)))
  getwd()
})

## ---------------------------------------------------------------------------
## Argument parsing
## ---------------------------------------------------------------------------

raw <- commandArgs(trailingOnly=TRUE)

if (length(raw) == 0 || any(raw %in% c("-h", "--help"))) {
  cat(usage_text())
  quit(save="no", status=if (length(raw) == 0) 1 else 0)
}

opts <- list(
  rds_files     = NULL,
  celltype_col  = NULL,
  out_dir       = NULL,
  label         = format(Sys.Date(), "%Y%m%d"),
  k1            = "10",
  k2            = NULL,
  k2_min        = "3",
  k2_max        = NULL,
  loo_mode      = "count",
  loo_n         = "1",
  agg_fn        = "median",
  condition_col = NULL,
  condition_map = NULL,
  seed          = "42",
  no_plots      = FALSE
)

flag_slot <- list(
  "--rds-files"     = "rds_files",
  "--celltype-col"  = "celltype_col",
  "--out-dir"       = "out_dir",
  "--label"         = "label",
  "--k1"            = "k1",
  "--k2"            = "k2",
  "--k2-min"        = "k2_min",
  "--k2-max"        = "k2_max",
  "--loo-mode"      = "loo_mode",
  "--loo-n"         = "loo_n",
  "--agg-fn"        = "agg_fn",
  "--condition-col" = "condition_col",
  "--condition-map" = "condition_map",
  "--seed"          = "seed"
)

has_flags <- any(grepl("^--", raw))

if (has_flags) {
  i <- 1L
  while (i <= length(raw)) {
    tok <- raw[i]
    if (tok == "--no-plots") { opts$no_plots <- TRUE; i <- i + 1L; next }
    if (!grepl("^--", tok))
      die("unexpected positional argument '", tok, "' in flag-style call (see --help)")
    if (grepl("=", tok, fixed=TRUE)) {
      key <- sub("=.*$",    "", tok)
      val <- sub("^[^=]*=", "", tok)
      i   <- i + 1L
    } else {
      key <- tok
      if (i + 1L > length(raw)) die("missing value for ", key)
      val <- raw[i + 1L]; i <- i + 2L
    }
    slot <- flag_slot[[key]]
    if (is.null(slot)) die("unknown option: ", key, " (see --help)")
    opts[[slot]] <- val
  }
} else {
  if (length(raw) != 11)
    die("positional form requires exactly 11 args; got ", length(raw),
        ". Run with --help.")
  opts$rds_files     <- raw[1]
  opts$celltype_col  <- raw[2]
  opts$out_dir       <- raw[3]
  opts$label         <- raw[4]
  opts$k1            <- raw[5]
  opts$k2_min        <- raw[6]
  opts$k2_max        <- raw[7]
  opts$loo_mode      <- raw[8]
  opts$loo_n         <- raw[9]
  opts$condition_map <- raw[10]
  opts$seed          <- raw[11]
}

## ---------------------------------------------------------------------------
## Normalize + validate
## ---------------------------------------------------------------------------

rds_paths <- parse_rds_files(opts$rds_files)
if (is.null(rds_paths)) die("--rds-files is required (see --help)")
missing_f <- rds_paths[!file.exists(rds_paths)]
if (length(missing_f))
  die("RDS file(s) not found:\n  ", paste(missing_f, collapse="\n  "))

if (is_unset(opts$out_dir)) {
  opts$out_dir <- getwd()
  message("run-neighborhood_analysis.R: no --out-dir given; using ", opts$out_dir)
}
dir.create(opts$out_dir, showWarnings=FALSE, recursive=TRUE)

sample_names <- tools::file_path_sans_ext(basename(rds_paths))
if (anyDuplicated(sample_names))
  sample_names <- paste0(sample_names, "_", seq_along(sample_names))
spe_paths_named <- setNames(rds_paths, sample_names)

label      <- trimws(opts$label %||% format(Sys.Date(), "%Y%m%d"))
k1         <- parse_int(opts$k1,    "--k1")
k2         <- if (is_unset(opts$k2)) NULL else parse_int(opts$k2, "--k2")
k2_min     <- parse_int(opts$k2_min, "--k2-min")
k2_max     <- if (is_unset(opts$k2_max)) NULL else parse_int(opts$k2_max, "--k2-max")
loo_n      <- parse_num(opts$loo_n, "--loo-n")
seed       <- parse_int(opts$seed,  "--seed")
loo_mode   <- opts$loo_mode
agg_fn     <- opts$agg_fn
cond_col   <- if (is_unset(opts$condition_col)) NULL else opts$condition_col
cond_map   <- parse_condition_map(opts$condition_map)
make_plots <- !isTRUE(opts$no_plots)

if (!loo_mode %in% c("count", "pct"))
  die("--loo-mode must be 'count' or 'pct' (got '", loo_mode, "')")
if (!agg_fn %in% c("median", "mean"))
  die("--agg-fn must be 'median' or 'mean' (got '", agg_fn, "')")

ct_cols_raw   <- parse_ct_cols(opts$celltype_col)
celltype_cols <- if (is.null(ct_cols_raw)) {
  NULL
} else if (is.null(names(ct_cols_raw))) {
  setNames(rep(ct_cols_raw, length(sample_names)), sample_names)
} else {
  ct_cols_raw
}

show <- function(lbl, value) {
  pretty <- if (is.null(value) || (length(value)==1 && is.na(value))) "<none>"
             else paste(value, collapse="; ")
  message(sprintf("  %-22s : %s", lbl, pretty))
}
message("== NeighborhoodR config ==")
show("rds_files",    paste(rds_paths, collapse="; "))
show("sample_names", sample_names)
show("out_dir",      opts$out_dir)
show("label",        label)
show("celltype_cols", if (is.null(celltype_cols)) "<auto>"
                      else paste(paste0(names(celltype_cols), "=",
                                        celltype_cols), collapse=", "))
show("k1",           k1)
show("k2",           k2 %||% "<sweep>")
show("k2_min/max",  paste0(k2_min, " – ", k2_max %||% "auto"))
show("loo_mode",     loo_mode)
show("loo_n",        loo_n)
show("agg_fn",       agg_fn)
show("condition_col", cond_col)
show("condition_map", if (is.null(cond_map)) "<none>"
                      else paste(paste0(names(cond_map), "=", cond_map),
                                 collapse=", "))
show("seed",         seed)
show("make_plots",   make_plots)
message("===========================")

## ---------------------------------------------------------------------------
## Source utilities
## ---------------------------------------------------------------------------

utils_file <- file.path(script_dir, "neighborhood_analysis-utils.R")
if (!file.exists(utils_file))
  die("neighborhood_analysis-utils.R not found under '", script_dir,
      "' (set NEIGHBORHOODR_DIR to override)")

Sys.setenv(NEIGHBORHOODR_DIR = script_dir)
source(utils_file)

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

run_neighborhood_analysis(
  spe_paths     = spe_paths_named,
  celltype_cols = celltype_cols,
  out_dir       = opts$out_dir,
  label         = label,
  k1            = k1,
  k2            = k2,
  k2_min        = k2_min,
  k2_max        = k2_max,
  loo_mode      = loo_mode,
  loo_n         = loo_n,
  agg_fn_name   = agg_fn,
  condition_col = cond_col,
  condition_map = cond_map,
  seed          = seed,
  make_plots    = make_plots,
  verbose       = TRUE
)
