#!/usr/bin/env Rscript
#
# run-pcf.R
# =========
# CLI entry point for the PCF (pair correlation function) pipeline — headless
# port of the PhenoSuite `pcf-v2` Shiny app.
#
# Usage:
#   Rscript run-pcf.R \
#       --vectra-files=DIR|a.csv;b.csv|list.txt \
#       --out-dir=DIR \
#       [--label=LABEL] \
#       [--celltypes=ct1,ct2,ct3] \
#       [--ref-celltype=ct1] \
#       [--radius=30] \
#       [--resolution=0.377] \
#       [--count-threshold=10] \
#       [--min-count=0] \
#       [--phenotype-col=Phenotype] \
#       [--no-plots]
#
# --vectra-files accepts:
#   * A directory — every *.csv beneath it is used (recursively), matching the
#     GUI's extract_data() glob.
#   * Semicolon-separated CSV paths:  a.csv;b.csv;c.csv
#   * A .txt file listing one CSV path per line.
#
# One CSV = one sample, as in the GUI.
#
# Run with --help for full details.

## ---------------------------------------------------------------------------
## Usage / help
## ---------------------------------------------------------------------------

usage_text <- function() {
"Usage:
  Rscript run-pcf.R --vectra-files=SPEC --out-dir=DIR [options]

Computes inhomogeneous pair correlation functions between cell types from
Vectra cell_seg_data-style CSVs, writes per-interaction PCF curves, a
pcf-builder-compatible AUC table, curve-grid plots and AUC violin plots.

Required:
  --vectra-files=SPEC   Directory of Vectra CSVs (searched recursively), OR
                        semicolon-separated CSV paths, OR a .txt file listing
                        one CSV path per line. One CSV = one sample.
  --out-dir=DIR         Root output directory (a label subdirectory is created).

Options:
  --label=STR           Output file prefix.                  [default: YYYYMMDD]
  --celltypes=LIST      Comma-separated cell types to analyse.
                        [default: phenotypes present in every input file]
  --ref-celltype=NAME   Reference cell type for the AUC violins.
                        [default: first analysed cell type]
  --radius=N            Maximum radius in microns.           [default: 30]
  --resolution=N        Instrument resolution, microns/pixel.[default: 0.377]
  --count-threshold=N   Minimum cells of a type before its interactions are
                        computed (must be >= 1).             [default: 10]
  --min-count=N         Drop interactions whose smaller cell count is <= this
                        when assembling curves.              [default: 0]
  --phenotype-col=NAME  Column holding the cell type labels. [default: Phenotype]
  --no-plots            Skip PDF plot generation.
  -h, --help            Show this help and exit.

Sentinels meaning 'not set': NULL, null, NA, none, None, '' (empty).

Environment variables:
  PCF_DIR              Directory of this script + pcf-utils.R (auto-detected).

Output (written to --out-dir/--label/):
  {label}_pcf_summary.csv            One row per sample x interaction (counts,
                                     PCFsum, normalization, skipped flag)
  {label}_pcf_curves.csv             Long-format curves: one row per
                                     (sample, interaction, radius step)
  {label}_ppc.rds                    Per-cell-type curve tables (GUI ppc.rds)
  {label}-PCF_AUCs.csv               normPCF series per partner cell type
                                     (input format for pcf-builder)
  {label}-PCF-plots.pdf              Curve grid, one panel per cell type
  {label}-PCF_AUC_violins.pdf        AUC violins (single-sample runs)
  {label}-PCF_AUC_violins-global.pdf AUC violins pooled over samples
  {label}-PCF_AUC_violins-bySample.pdf   AUC violins split by sample
  individual-samples/*.pdf           Per-sample AUC violins (multi-sample runs)
  {label}_provenance.json            Full run provenance

Prerequisites:
  conda env create -f environment.yml && conda activate runpcf
"
}

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

die <- function(...) {
  message("run-pcf.R: error: ", ...)
  quit(save = "no", status = 1)
}

# Defined here as well as in pcf-utils.R: the config banner below prints
# before the utils are sourced (base R only gained %||% in 4.4).
`%||%` <- function(a, b) if (!is.null(a) && length(a) && !all(is.na(a))) a else b

is_unset <- function(x) {
  if (is.null(x) || length(x) == 0) return(TRUE)
  if (length(x) > 1)                return(FALSE)
  if (is.na(x))                     return(TRUE)
  trimws(as.character(x)) %in%
    c("", "NULL", "null", "NA", "na", "none", "None", "NONE")
}

parse_list <- function(x, sep = ",") {
  if (is_unset(x)) return(NULL)
  parts <- trimws(unlist(strsplit(as.character(x), sep, fixed = TRUE)))
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

# Directory | ;-separated paths | .txt list — mirrors the GUI's file chooser,
# which accepted a whole directory of Vectra exports.
parse_vectra_files <- function(x) {
  if (is_unset(x)) return(NULL)
  x <- trimws(as.character(x))
  if (dir.exists(x)) {
    files <- list.files(x, pattern = "\\.csv$", recursive = TRUE,
                        full.names = TRUE, ignore.case = TRUE)
    return(sort(files))
  }
  if (grepl("\\.txt$", x, ignore.case = TRUE) && file.exists(x)) {
    paths <- trimws(readLines(x, warn = FALSE))
    return(paths[nzchar(paths)])
  }
  parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  parts[nzchar(parts)]
}

## ---------------------------------------------------------------------------
## Script-directory discovery
## ---------------------------------------------------------------------------

script_dir <- local({
  env_dir <- Sys.getenv("PCF_DIR", unset = "")
  if (nzchar(env_dir)) {
    if (!dir.exists(env_dir)) die("PCF_DIR does not exist: ", env_dir)
    return(normalizePath(env_dir))
  }
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- sub("^--file=", "", grep("^--file=", args, value = TRUE))
  if (length(file_arg) == 1 && file.exists(file_arg))
    return(dirname(normalizePath(file_arg)))
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
  vectra_files    = NULL,
  out_dir         = NULL,
  label           = format(Sys.Date(), "%Y%m%d"),
  celltypes       = NULL,
  ref_celltype    = NULL,
  radius          = "30",
  resolution      = "0.377",
  count_threshold = "10",
  min_count       = "0",
  phenotype_col   = "Phenotype",
  no_plots        = FALSE
)

flag_slot <- list(
  "--vectra-files"    = "vectra_files",
  "--out-dir"         = "out_dir",
  "--label"           = "label",
  "--celltypes"       = "celltypes",
  "--ref-celltype"    = "ref_celltype",
  "--radius"          = "radius",
  "--resolution"      = "resolution",
  "--count-threshold" = "count_threshold",
  "--min-count"       = "min_count",
  "--phenotype-col"   = "phenotype_col"
)

i <- 1L
while (i <= length(raw)) {
  tok <- raw[i]
  if (tok == "--no-plots") { opts$no_plots <- TRUE; i <- i + 1L; next }
  if (!grepl("^--", tok))
    die("unexpected positional argument '", tok, "' (see --help)")
  if (grepl("=", tok, fixed = TRUE)) {
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

## ---------------------------------------------------------------------------
## Normalize + validate
## ---------------------------------------------------------------------------

vectra_files <- parse_vectra_files(opts$vectra_files)
if (is.null(vectra_files) || length(vectra_files) == 0)
  die("--vectra-files is required and must resolve to at least one CSV ",
      "(see --help)")
missing_f <- vectra_files[!file.exists(vectra_files)]
if (length(missing_f))
  die("Vectra file(s) not found:\n  ", paste(missing_f, collapse = "\n  "))

if (is_unset(opts$out_dir)) {
  opts$out_dir <- getwd()
  message("run-pcf.R: no --out-dir given; using ", opts$out_dir)
}
dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)

label           <- if (is_unset(opts$label)) format(Sys.Date(), "%Y%m%d") else
                     trimws(opts$label)
cell_types      <- parse_list(opts$celltypes)
ref_celltype    <- if (is_unset(opts$ref_celltype)) NULL else trimws(opts$ref_celltype)
radius          <- parse_num(opts$radius,     "--radius")
resolution      <- parse_num(opts$resolution, "--resolution")
count_threshold <- parse_int(opts$count_threshold, "--count-threshold")
min_count       <- parse_num(opts$min_count,  "--min-count")
phenotype_col   <- if (is_unset(opts$phenotype_col)) "Phenotype" else
                     trimws(opts$phenotype_col)
make_plots      <- !isTRUE(opts$no_plots)

if (radius <= 0)          die("--radius must be > 0 (got ", radius, ")")
if (resolution <= 0)      die("--resolution must be > 0 (got ", resolution, ")")
if (count_threshold < 1L) die("--count-threshold must be >= 1 (got ",
                              count_threshold, ")")

show <- function(lbl, value) {
  pretty <- if (is.null(value) || (length(value) == 1 && is.na(value))) "<none>"
            else paste(value, collapse = "; ")
  message(sprintf("  %-18s : %s", lbl, pretty))
}
message("== PCF config ==")
show("vectra_files",    paste(basename(vectra_files), collapse = "; "))
show("n_samples",       length(vectra_files))
show("out_dir",         opts$out_dir)
show("label",           label)
show("celltypes",       cell_types %||% "<shared across inputs>")
show("ref_celltype",    ref_celltype %||% "<first celltype>")
show("radius (um)",     radius)
show("resolution",      resolution)
show("count_threshold", count_threshold)
show("min_count",       min_count)
show("phenotype_col",   phenotype_col)
show("make_plots",      make_plots)
message("================")

## ---------------------------------------------------------------------------
## Source utilities
## ---------------------------------------------------------------------------

utils_file <- file.path(script_dir, "pcf-utils.R")
if (!file.exists(utils_file))
  die("pcf-utils.R not found under '", script_dir, "' (set PCF_DIR to override)")

Sys.setenv(PCF_DIR = script_dir)
source(utils_file)

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

run_pcf(
  vectra_files    = vectra_files,
  out_dir         = opts$out_dir,
  label           = label,
  cell_types      = cell_types,
  ref_celltype    = ref_celltype,
  radius          = radius,
  resolution      = resolution,
  count_threshold = count_threshold,
  min_count       = min_count,
  phenotype_col   = phenotype_col,
  make_plots      = make_plots,
  verbose         = TRUE
)
