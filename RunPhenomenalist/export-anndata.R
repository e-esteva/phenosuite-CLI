#!/usr/bin/env Rscript
#
# export-anndata.R
# =================
# Standalone converter: reads an existing spe.rds (SpatialExperiment, as
# written by RunPhenomenalist / run-phenomenalist.R) and writes it out as an
# AnnData .h5ad file, without re-running the clustering pipeline. Useful for
# samples processed before --export-anndata existed on run-phenomenalist.R,
# or for one-off conversions.
#
# Usage:
#   Rscript export-anndata.R --spe-rds=PATH [--out=PATH] [--assay=NAME]
#
# Run with --help for details.

## ---------------------------------------------------------------------------
## Usage / help
## ---------------------------------------------------------------------------

usage_text <- function() {
"Usage:
  Rscript export-anndata.R --spe-rds=PATH [options]

Converts an existing spe.rds (SpatialExperiment) to an AnnData .h5ad file.
Cell coordinates (spatialCoords) are carried over into adata.obsm['spatial'].

Required:
  --spe-rds=PATH   Path to a spe.rds file written by RunPhenomenalist.

Options:
  --out=PATH       Output .h5ad path.            [default: same dir/name as --spe-rds, .h5ad extension]
  --assay=NAME     Assay to store as AnnData X.   [default: counts, falls back to the first assay]
  -h, --help       Show this help and exit.

Environment variables:
  PHENOMENALIST_DIR   Directory holding phenomenalist-utils.R. Defaults to this
                       script's own directory.

Prerequisites:
  Requires the Bioconductor 'zellkonverter' package:
      Rscript -e 'BiocManager::install(\"zellkonverter\")'
"
}

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

die <- function(...) {
  message("export-anndata.R: error: ", ...)
  quit(save = "no", status = 1)
}

raw <- commandArgs(trailingOnly = TRUE)

if (length(raw) == 0 || any(raw %in% c("-h", "--help"))) {
  cat(usage_text())
  quit(save = "no", status = if (length(raw) == 0) 1 else 0)
}

## ---------------------------------------------------------------------------
## Script-directory discovery (for sourcing phenomenalist-utils.R portably)
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

opts <- list(spe_rds = NULL, out = NULL, assay = "counts")

flag_slot <- list(
  "--spe-rds" = "spe_rds",
  "--out"     = "out",
  "--assay"   = "assay"
)

i <- 1L
while (i <= length(raw)) {
  tok <- raw[i]
  if (!grepl("^--", tok)) die("unexpected argument: '", tok, "' (see --help)")
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

## ---------------------------------------------------------------------------
## Normalize + validate
## ---------------------------------------------------------------------------

if (is.null(opts$spe_rds) || !nzchar(opts$spe_rds)) die("--spe-rds is required")
if (!file.exists(opts$spe_rds)) die("spe.rds not found: ", opts$spe_rds)

out_path <- opts$out
if (is.null(out_path) || !nzchar(out_path)) {
  out_path <- sub("\\.rds$", ".h5ad", opts$spe_rds)
  if (identical(out_path, opts$spe_rds)) out_path <- paste0(opts$spe_rds, ".h5ad")
}
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

message("== export-anndata.R config ==")
message("  spe_rds : ", opts$spe_rds)
message("  out     : ", out_path)
message("  assay   : ", opts$assay)
message("==============================")

## ---------------------------------------------------------------------------
## Load dependencies + phenomenalist-utils.R (for export_anndata.mod)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  suppressWarnings({
    library(glue)
    library(SpatialExperiment)
  })
})

utils_file <- file.path(script_dir, "phenomenalist-utils.R")
if (!file.exists(utils_file)) {
  die("phenomenalist-utils.R not found under '", script_dir, "' (set PHENOMENALIST_DIR to override)")
}
source(utils_file)

## ---------------------------------------------------------------------------
## Convert
## ---------------------------------------------------------------------------

message("reading ", opts$spe_rds)
spe <- readRDS(opts$spe_rds)

export_anndata.mod(
  spe,
  out_dir  = dirname(out_path),
  filename = basename(out_path),
  X_name   = opts$assay
)

message("wrote ", out_path)
