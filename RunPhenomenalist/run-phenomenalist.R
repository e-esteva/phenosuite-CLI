#!/usr/bin/env Rscript
#
# run-phenomenalist.R
# ===================
# CLI entry point for the RunPhenomenalist pipeline.
#
# Supports two invocation styles:
#
#   (1) Named flags (preferred):
#       Rscript run-phenomenalist.R \
#           --segmentation-file=PATH \
#           --out-dir=DIR \
#           [--failed-markers=LIST] \
#           [--nuclear-markers=LIST] \
#           [--classifier-label=LIST] \
#           [--clustering-res=SPEC] \
#           [--max-cells=N] \
#           [--phenotyping-template=PATH] \
#           [--skip-cols=REGEX]
#
#   (2) Positional (legacy, for run-phenomenalist.s back-compat):
#       Rscript run-phenomenalist.R \
#           <segmentation_file> <failed_markers> <nuclear_markers> \
#           <out_dir> <clustering_res> <classifier_label> \
#           <max_cells> <phenotyping_template>
#
# Segmentation format (HALO / Mesmer / QuPath) is auto-detected from the
# column headers; QuPath centroid variants ("Centroid X µm", "Centroid X px",
# "Centroid X") are all resolved to x/y automatically. Pass --skip-cols to
# override the column-filter regex for custom schemas or when auto-detection
# picks wrong.
#
# Run with --help for details.

## ---------------------------------------------------------------------------
## Usage / help
## ---------------------------------------------------------------------------

usage_text <- function() {
"Usage:
  Rscript run-phenomenalist.R --segmentation-file=PATH --out-dir=DIR [options]
  Rscript run-phenomenalist.R <seg> <failed> <nuclear> <out> <res> <label> <max> <template>

Runs the RunPhenomenalist phenotyping + clustering pipeline on one segmentation file.
Segmentation format (HALO / Mesmer / QuPath) is auto-detected from the column headers.
QuPath centroid variants ('Centroid X µm', 'Centroid X px', 'Centroid X') all resolve to x/y.

Required:
  --segmentation-file=PATH    Per-cell segmentation CSV / TSV / TSV.gz.
  --out-dir=DIR               Output directory (created if missing).

Options:
  --failed-markers=LIST       Comma-separated markers to drop.        [default: none]
  --nuclear-markers=LIST      Comma-separated nuclear markers.        [default: none]
  --classifier-label=LIST     Comma-separated classifier labels.      [default: none]
  --clustering-res=SPEC       Clustering resolutions. Accepts:
                                * comma list  : '1,2,3'  -> c(1,2,3)
                                * range       : '5:7'    -> c(5,6,7)
                              [default: 1,2]
  --max-cells=N               Cell subsample threshold.               [default: 100000]
  --phenotyping-template=PATH Manual-gating template CSV.             [default: none]
  --skip-cols=REGEX           Override the column-filter regex used to exclude
                              non-marker columns. Use this when the segmentation
                              schema is non-standard or auto-detection is wrong.
                                                                      [default: auto]
  -h, --help                  Show this help and exit.

Sentinels that all mean 'not set' for any option value:
  NULL, null, NA, none, None, '' (empty), 0.

Environment variables:
  PHENOMENALIST_DIR           Directory holding RunPhenomenalist.R and phenomenalist-utils.R.
                              Defaults to the directory containing this script.

Prerequisites:
  The 'phenomenalist' R package must be installed (loaded via library()). Install from
  GitHub with:
      Rscript -e 'remotes::install_github(\"igordot/phenomenalist\")'
  Or provision the shipped conda env:
      conda env create -f environment.yml && conda activate runphenomenalist
      Rscript -e 'remotes::install_github(\"igordot/phenomenalist\")'
"
}

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

die <- function(...) {
  message("run-phenomenalist.R: error: ", ...)
  quit(save = "no", status = 1)
}

is_unset <- function(x) {
  if (is.null(x)) return(TRUE)
  if (length(x) == 0) return(TRUE)
  if (length(x) > 1) return(FALSE)
  if (is.na(x)) return(TRUE)
  trimws(as.character(x)) %in% c("", "0", "NULL", "null", "NA", "none", "None", "NONE")
}

parse_list <- function(x) {
  if (is_unset(x)) return(NULL)
  parts <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) NULL else parts
}

parse_res <- function(x) {
  if (is_unset(x)) return(c(1, 2))
  x <- trimws(as.character(x))
  if (grepl(":", x, fixed = TRUE)) {
    bounds <- suppressWarnings(as.numeric(strsplit(x, ":", fixed = TRUE)[[1]]))
    if (length(bounds) != 2 || any(is.na(bounds))) {
      die("--clustering-res range must be 'low:high' (got '", x, "')")
    }
    return(seq(bounds[1], bounds[2]))
  }
  vals <- suppressWarnings(as.numeric(unlist(strsplit(x, ",", fixed = TRUE))))
  if (any(is.na(vals))) die("--clustering-res list must be numeric (got '", x, "')")
  vals
}

parse_numeric <- function(x, flag) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) die(flag, " must be numeric (got '", x, "')")
  v
}

## ---------------------------------------------------------------------------
## Script-directory discovery (for sourcing sibling files portably)
## ---------------------------------------------------------------------------

script_dir <- local({
  env_dir <- Sys.getenv("PHENOMENALIST_DIR", unset = "")
  if (nzchar(env_dir)) {
    if (!dir.exists(env_dir)) die("PHENOMENALIST_DIR is set but does not exist: ", env_dir)
    return(normalizePath(env_dir))
  }
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file_arg) == 1 && file.exists(file_arg)) {
    return(dirname(normalizePath(file_arg)))
  }
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

opts <- list(
  segmentation_file    = NULL,
  failed_markers       = NULL,
  nuclear_markers      = NULL,
  out_dir              = NULL,
  clustering_res       = "1,2",
  classifier_label     = NULL,
  max_cells            = "100000",
  phenotyping_template = NULL,
  skip_cols            = NULL
)

flag_slot <- list(
  "--segmentation-file"    = "segmentation_file",
  "--failed-markers"       = "failed_markers",
  "--nuclear-markers"      = "nuclear_markers",
  "--out-dir"              = "out_dir",
  "--clustering-res"       = "clustering_res",
  "--classifier-label"     = "classifier_label",
  "--max-cells"            = "max_cells",
  "--phenotyping-template" = "phenotyping_template",
  "--skip-cols"            = "skip_cols"
)

has_flags <- any(grepl("^--", raw))

if (has_flags) {
  i <- 1L
  while (i <= length(raw)) {
    tok <- raw[i]
    if (!grepl("^--", tok)) {
      die("unexpected positional argument in flag-style call: '", tok, "' (see --help)")
    }
    if (grepl("=", tok, fixed = TRUE)) {
      key <- sub("=.*$", "", tok)
      val <- sub("^[^=]*=", "", tok)
      i <- i + 1L
    } else {
      key <- tok
      if (i + 1L > length(raw)) die("missing value for ", key)
      val <- raw[i + 1L]
      i <- i + 2L
    }
    slot <- flag_slot[[key]]
    if (is.null(slot)) die("unknown option: ", key, " (see --help)")
    opts[[slot]] <- val
  }
} else {
  if (length(raw) != 8) {
    die("positional form requires exactly 8 args; got ", length(raw), ". Run with --help.")
  }
  opts$segmentation_file    <- raw[1]
  opts$failed_markers       <- raw[2]
  opts$nuclear_markers      <- raw[3]
  opts$out_dir              <- raw[4]
  opts$clustering_res       <- raw[5]
  opts$classifier_label     <- raw[6]
  opts$max_cells            <- raw[7]
  opts$phenotyping_template <- raw[8]
}

## ---------------------------------------------------------------------------
## Normalize + validate
## ---------------------------------------------------------------------------

if (is_unset(opts$segmentation_file)) die("--segmentation-file is required")
if (!file.exists(opts$segmentation_file)) {
  die("segmentation file not found: ", opts$segmentation_file)
}

if (is_unset(opts$out_dir)) {
  opts$out_dir <- getwd()
  message("run-phenomenalist.R: no --out-dir given; using ", opts$out_dir)
}
dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)

max_cells            <- parse_numeric(opts$max_cells, "--max-cells")
clustering_res       <- parse_res(opts$clustering_res)
failed_markers       <- parse_list(opts$failed_markers)
nuclear_markers      <- parse_list(opts$nuclear_markers)
classifier_label     <- parse_list(opts$classifier_label)
phenotyping_template <- if (is_unset(opts$phenotyping_template)) NULL else opts$phenotyping_template
skip_cols            <- if (is_unset(opts$skip_cols)) NULL else as.character(opts$skip_cols)

if (!is.null(phenotyping_template) && !file.exists(phenotyping_template)) {
  die("phenotyping template not found: ", phenotyping_template)
}

show <- function(label, value) {
  pretty <- if (is.null(value)) "<none>" else paste(value, collapse = ",")
  message(sprintf("  %-22s : %s", label, pretty))
}
message("== RunPhenomenalist config ==")
show("segmentation_file",    opts$segmentation_file)
show("out_dir",              opts$out_dir)
show("max_cells",            max_cells)
show("clustering_res",       clustering_res)
show("failed_markers",       failed_markers)
show("nuclear_markers",      nuclear_markers)
show("classifier_label",     classifier_label)
show("phenotyping_template", phenotyping_template)
show("skip_cols",            skip_cols)
message("=============================")

## ---------------------------------------------------------------------------
## Load pipeline dependencies (deferred so --help / arg errors stay clean)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  suppressWarnings({
    library(glue)
    library(MatrixGenerics)
    library(SpatialExperiment)
  })
})

## ---------------------------------------------------------------------------
## Source pipeline code
## ---------------------------------------------------------------------------

pipeline_file <- file.path(script_dir, "RunPhenomenalist.R")
if (!file.exists(pipeline_file)) {
  die("RunPhenomenalist.R not found under '", script_dir, "' (set PHENOMENALIST_DIR to override)")
}
utils_file <- file.path(script_dir, "phenomenalist-utils.R")
if (!file.exists(utils_file)) {
  die("phenomenalist-utils.R not found under '", script_dir, "' (set PHENOMENALIST_DIR to override)")
}

# Export for RunPhenomenalist.R so its own source() line can find the utils.
Sys.setenv(PHENOMENALIST_DIR = script_dir)

source(utils_file)
source(pipeline_file)

## ---------------------------------------------------------------------------
## Load the public phenomenalist package (https://github.com/igordot/phenomenalist)
## ---------------------------------------------------------------------------

ok <- tryCatch({
  suppressPackageStartupMessages(suppressWarnings(library(phenomenalist)))
  TRUE
}, error = function(e) FALSE)
if (!ok) {
  die("phenomenalist package is not installed. Install from GitHub with:\n",
      "  Rscript -e 'remotes::install_github(\"igordot/phenomenalist\")'\n",
      "Or provision the shipped conda env: conda env create -f environment.yml")
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

RunPhenomenalist(
  segmentation_file    = opts$segmentation_file,
  failed.markers       = failed_markers,
  nuclear.markers      = nuclear_markers,
  out_dir              = opts$out_dir,
  clustering_res       = clustering_res,
  classifier_label     = classifier_label,
  max.cells            = max_cells,
  min.cells            = 10,
  phenotyping_template = phenotyping_template,
  skip_cols            = skip_cols
)
