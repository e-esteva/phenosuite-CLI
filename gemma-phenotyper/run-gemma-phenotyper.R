#!/usr/bin/env Rscript
#
# run-gemma-phenotyper.R
# ======================
# CLI entry point for the Gemma mIF phenotyper pipeline.
#
#   Rscript run-gemma-phenotyper.R \
#       --spe-file=PATH \
#       --out-dir=DIR \
#       [--backend=hf|ollama] [--prompt-style=zero-shot|ontology] \
#       [--cluster-col=NAME] [--markers=LIST] [--mode=cluster|cell] \
#       [--ollama-host=URL] [--ollama-model=TAG] [--temperature=F] [--harmonize=BOOL] \
#       [--model-dir=PATH] [--base-model=ID] [--load-in-8bit=BOOL] [--ontology=PATH] \
#       [--tissue=STR] [--label=STR] [--export-loom=BOOL]
#
# Two backends:
#   hf (default) — HuggingFace-format weights run through infer.py
#                  (transformers/PEFT). Works offline on an HPC node once the
#                  weights are on disk. Requires --model-dir.
#   ollama       — HTTP call to a running Ollama server. Requires a reachable
#                  server, so usually only viable off-cluster.
#
# Prompt style is independent of the backend: --prompt-style=zero-shot uses the
# tissue-aware "most distinctive markers" prompt (for models not fine-tuned on
# your ontology); --prompt-style=ontology sends every marker plus an explicit
# ontology table (for checkpoints fine-tuned against it).
#
# Run with --help for details.

## ---------------------------------------------------------------------------
## Usage / help
## ---------------------------------------------------------------------------

usage_text <- function() {
"Usage:
  Rscript run-gemma-phenotyper.R --spe-file=PATH --out-dir=DIR --model-dir=PATH [options]
  Rscript run-gemma-phenotyper.R --spe-file=PATH --out-dir=DIR --backend=ollama --ollama-model=TAG [options]

Annotates clusters (or cells) in a SpatialExperiment / SingleCellExperiment with
LLM-assigned cell types, using a Gemma model via Ollama or a fine-tuned checkpoint.

Required:
  --spe-file=PATH             SpatialExperiment / SingleCellExperiment .rds.
  --out-dir=DIR               Output directory (created if missing).

Backend selection:
  --backend=NAME              'hf' (HuggingFace weights via infer.py — works
                              offline, the HPC option) or 'ollama' (HTTP call to
                              a running Ollama server).                [default: hf]
  --prompt-style=NAME         'zero-shot' (tissue + top/bottom 10% most
                              distinctive markers; for models NOT fine-tuned on
                              your ontology) or 'ontology' (all markers plus an
                              explicit ontology table; for fine-tuned
                              checkpoints).                            [default: zero-shot]

HF backend (--backend=hf):
  --model-dir=PATH            Merged checkpoint or LoRA adapter directory.
                              Type is auto-detected from config.json /
                              adapter_config.json.                     [required]
  --base-model=ID             Base model (HF id or local path); only used when
                              --model-dir is a LoRA adapter.
                                                                       [default: google/gemma-3-1b-it]
  --load-in-8bit=BOOL         Load in 8-bit (needs bitsandbytes).      [default: false]
  --python=PATH               Python interpreter that has transformers/PEFT.
                              [default: $GEMMA_PYTHON or python3]

Ollama backend (--backend=ollama):
  --ollama-model=TAG          Model tag, e.g. 'gemma3:12b-it-qat'.     [required]
  --ollama-host=URL           Ollama server URL.
                              [default: $OLLAMA_HOST or http://localhost:11434]
  --temperature=F             Sampling temperature, 0-1. 0 is greedy decoding:
                              identical input always yields the same label,
                              which is what a reproducible batch run needs.
                              Raise it only to sample one cluster repeatedly and
                              use the agreement rate as an empirical confidence
                              (see vote_unclassified_ollama).            [default: 0]
  --timeout=N                 Per-request timeout in seconds. Large models on
                              CPU-only nodes can take minutes per call, and a
                              blown timeout is recorded as a 'request_failure'
                              prediction rather than aborting the run.
                                                                       [default: 600]
  --marker-glossary=PATH      Glossary file ('MARKER = gloss' per line) telling the
                              model what your panel's non-obvious markers mean.
                              Only entries for a cluster's elevated markers are
                              injected. Strongly recommended: without it a
                              zero-shot model guesses at transcription factors and
                              alias names, confidently and wrongly.
                              See marker-glossary-example.txt.       [default: none]
  --min-z=N                   A marker counts as 'elevated' only if it exceeds
                              this value. Clusters with nothing above it are
                              sent as 'no markers are elevated' and answered
                              Unclassified rather than guessed. Only used with
                              --prompt-style=zero-shot.                [default: 0.5]
  --vote-samples=N            Re-sample each Unclassified cluster N times and use
                              the agreement rate as an empirical confidence (the
                              model's self-reported confidence is not
                              informative). Runs in-process on both backends, so
                              the checkpoint is not reloaded. 0 disables. [default: 0]
  --vote-temperature=F        Temperature for the vote. Must be > 0 or every
                              sample is identical.                     [default: 0.7]
  --vote-min-agreement=F      Fraction of votes one label needs before it replaces
                              Unclassified.                            [default: 0.6]
  --vote-override-floor=BOOL  Also let voting override clusters that had NO
                              elevated markers. Off by default: a majority there
                              is the model guessing from nothing.      [default: false]
  --retry=N                   Extra passes over clusters that produced no usable
                              label (request or parse failure). Runs in-process
                              on both backends, so the HF checkpoint is not
                              reloaded. 0 disables.                    [default: 1]
  --harmonize=BOOL            Run the label-harmonisation pass that collapses
                              near-duplicate labels ('Treg cell' / 'T regulatory'
                              -> one canonical form). Works on both backends: an
                              extra HTTP call for ollama, or an extra generation
                              inside the same infer.py process for hf (so the
                              checkpoint is loaded only once).         [default: true]

Common options:
  --ontology=PATH             Text file holding the marker -> cell-type ontology.
                              Only read when --prompt-style=ontology.
                                                                       [default: built-in]
  --mode=NAME                 'cluster' (one prompt per cluster, fast) or
                              'cell' (one prompt per cell, slow).      [default: cluster]
  --cluster-col=NAME          colData column holding cluster IDs.      [required in cluster mode]
  --markers=LIST              Comma-separated features to include.     [default: all rownames]
  --tissue=STR                Tissue context, e.g. 'human lymph node'. Added to
                              the prompt and used as the Vectra tissue fallback.
                                                                       [default: none]
  --label=STR                 Output filename prefix.                  [default: .rds basename]
  --export-loom=BOOL          Also write a .loom alongside the .rds.   [default: false]
  -h, --help                  Show this help and exit.

Sentinels that all mean 'not set' for any option value:
  NULL, null, NA, none, None, '' (empty).

Environment variables:
  GEMMA_DIR                   Directory holding RunGemmaPhenotyper.R, gemma-utils.R,
                              and infer.py. Defaults to this script's directory.
  OLLAMA_HOST                 Default Ollama URL when --ollama-host is omitted.
  GEMMA_PYTHON                Default Python interpreter for the local backend.

Outputs (prefixed with --label):
  *_spe_gemma_annotated.rds       annotated object (gemma_cell_type/gemma_confidence)
  *_gemma_predictions.csv         per-cluster/per-cell raw model output
  *_gemma_colData.csv             full colData dump
  *_gemma_celltype-summary.csv    cell counts per assigned type
  *_harmonisation-map.csv         raw -> canonical label mapping (ollama + --harmonize)
  vectra_gemma_*.csv              Vectra-format export (SpatialExperiment only)
  *_gemma_annotated.loom          optional, with --export-loom=true

Prerequisites:
  conda env create -f environment.yml && conda activate gemma-phenotyper
  Ollama backend additionally needs a reachable server:  ollama serve
                                                         ollama pull gemma3:12b-it-qat
"
}

## ---------------------------------------------------------------------------
## Small helpers
## ---------------------------------------------------------------------------

die <- function(...) {
  message("run-gemma-phenotyper.R: error: ", ...)
  quit(save = "no", status = 1)
}

is_unset <- function(x) {
  if (is.null(x)) return(TRUE)
  if (length(x) == 0) return(TRUE)
  if (length(x) > 1) return(FALSE)
  if (is.na(x)) return(TRUE)
  trimws(as.character(x)) %in% c("", "NULL", "null", "NA", "none", "None", "NONE")
}

parse_list <- function(x) {
  if (is_unset(x)) return(NULL)
  parts <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE)))
  parts <- parts[nzchar(parts)]
  if (length(parts) == 0) NULL else parts
}

parse_numeric <- function(x, flag) {
  v <- suppressWarnings(as.numeric(x))
  if (is.na(v)) die(flag, " must be numeric (got '", x, "')")
  v
}

parse_bool <- function(x, flag, default = FALSE) {
  if (is_unset(x)) return(default)
  v <- tolower(trimws(as.character(x)))
  if (v %in% c("true", "t", "yes", "y", "1")) return(TRUE)
  if (v %in% c("false", "f", "no", "n", "0")) return(FALSE)
  die(flag, " must be true/false (got '", x, "')")
}

parse_choice <- function(x, flag, choices, default) {
  if (is_unset(x)) return(default)
  v <- tolower(trimws(as.character(x)))
  if (!v %in% choices)
    die(flag, " must be one of ", paste(choices, collapse = "/"), " (got '", x, "')")
  v
}

## ---------------------------------------------------------------------------
## Script-directory discovery (for sourcing sibling files portably)
## ---------------------------------------------------------------------------

script_dir <- local({
  env_dir <- Sys.getenv("GEMMA_DIR", unset = "")
  if (nzchar(env_dir)) {
    if (!dir.exists(env_dir)) die("GEMMA_DIR is set but does not exist: ", env_dir)
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
  spe_file     = NULL,
  out_dir      = NULL,
  backend      = "hf",
  prompt_style = "zero-shot",
  mode         = "cluster",
  cluster_col  = NULL,
  markers      = NULL,
  tissue       = NULL,
  label        = NULL,
  ollama_host  = NULL,
  ollama_model = NULL,
  temperature  = "0",
  timeout      = "600",
  retry        = "1",
  min_z        = "0.5",
  marker_glossary = NULL,
  vote_samples = "0",
  vote_temperature = "0.7",
  vote_min_agreement = "0.6",
  vote_override_floor = "false",
  harmonize    = "true",
  model_dir    = NULL,
  base_model   = "google/gemma-3-1b-it",
  load_in_8bit = "false",
  ontology     = NULL,
  python       = NULL,
  export_loom  = "false"
)

flag_slot <- list(
  "--spe-file"     = "spe_file",
  "--out-dir"      = "out_dir",
  "--backend"      = "backend",
  "--prompt-style" = "prompt_style",
  "--mode"         = "mode",
  "--cluster-col"  = "cluster_col",
  "--markers"      = "markers",
  "--tissue"       = "tissue",
  "--label"        = "label",
  "--ollama-host"  = "ollama_host",
  "--ollama-model" = "ollama_model",
  "--temperature"  = "temperature",
  "--timeout"      = "timeout",
  "--retry"        = "retry",
  "--min-z"        = "min_z",
  "--marker-glossary" = "marker_glossary",
  "--vote-samples" = "vote_samples",
  "--vote-temperature" = "vote_temperature",
  "--vote-min-agreement" = "vote_min_agreement",
  "--vote-override-floor" = "vote_override_floor",
  "--harmonize"    = "harmonize",
  "--model-dir"    = "model_dir",
  "--base-model"   = "base_model",
  "--load-in-8bit" = "load_in_8bit",
  "--ontology"     = "ontology",
  "--python"       = "python",
  "--export-loom"  = "export_loom"
)

i <- 1L
while (i <= length(raw)) {
  tok <- raw[i]
  if (!grepl("^--", tok)) die("unexpected positional argument: '", tok, "' (see --help)")
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

if (is_unset(opts$spe_file)) die("--spe-file is required")
if (!file.exists(opts$spe_file)) die("spe file not found: ", opts$spe_file)

if (is_unset(opts$out_dir)) {
  opts$out_dir <- getwd()
  message("run-gemma-phenotyper.R: no --out-dir given; using ", opts$out_dir)
}
dir.create(opts$out_dir, showWarnings = FALSE, recursive = TRUE)

# 'local' accepted as an alias for 'hf' — the flag was named that in the first cut.
if (!is_unset(opts$backend) && tolower(trimws(opts$backend)) == "local") opts$backend <- "hf"
backend      <- parse_choice(opts$backend, "--backend", c("hf", "ollama"), "hf")
prompt_style <- parse_choice(opts$prompt_style, "--prompt-style",
                             c("zero-shot", "ontology"), "zero-shot")
mode         <- parse_choice(opts$mode,    "--mode",    c("cluster", "cell"), "cluster")
markers      <- parse_list(opts$markers)
cluster_col  <- if (is_unset(opts$cluster_col)) NULL else trimws(opts$cluster_col)
tissue       <- if (is_unset(opts$tissue))      NULL else trimws(opts$tissue)
label        <- if (is_unset(opts$label))       NULL else trimws(opts$label)
harmonize    <- parse_bool(opts$harmonize,    "--harmonize",    TRUE)
load_in_8bit <- parse_bool(opts$load_in_8bit, "--load-in-8bit", FALSE)
export_loom  <- parse_bool(opts$export_loom,  "--export-loom",  FALSE)
temperature  <- parse_numeric(opts$temperature, "--temperature")
if (temperature < 0 || temperature > 1) die("--temperature must be between 0 and 1")
timeout      <- parse_numeric(opts$timeout, "--timeout")
retry        <- parse_numeric(opts$retry, "--retry")
min_z        <- parse_numeric(opts$min_z, "--min-z")
marker_glossary <- if (is_unset(opts$marker_glossary)) NULL else opts$marker_glossary
if (!is.null(marker_glossary) && !file.exists(marker_glossary))
  die("marker glossary not found: ", marker_glossary)
vote_samples <- parse_numeric(opts$vote_samples, "--vote-samples")
vote_temperature   <- parse_numeric(opts$vote_temperature, "--vote-temperature")
vote_min_agreement <- parse_numeric(opts$vote_min_agreement, "--vote-min-agreement")
vote_override_floor <- parse_bool(opts$vote_override_floor, "--vote-override-floor")
if (vote_samples < 0) die("--vote-samples must be >= 0")
if (vote_samples > 1 && vote_temperature <= 0)
  die("--vote-temperature must be > 0 when voting (0 returns identical samples)")
if (retry < 0) die("--retry must be >= 0")
if (timeout <= 0) die("--timeout must be a positive number of seconds")

if (mode == "cluster" && is.null(cluster_col))
  die("--cluster-col is required when --mode=cluster")

ollama_host  <- if (is_unset(opts$ollama_host))
  Sys.getenv("OLLAMA_HOST", "http://localhost:11434") else trimws(opts$ollama_host)
ollama_model <- if (is_unset(opts$ollama_model)) NULL else trimws(opts$ollama_model)
model_dir    <- if (is_unset(opts$model_dir))    NULL else trimws(opts$model_dir)
python_bin   <- if (is_unset(opts$python))
  Sys.getenv("GEMMA_PYTHON", "python3") else trimws(opts$python)

if (backend == "ollama") {
  if (is.null(ollama_model)) die("--ollama-model is required with --backend=ollama")
} else {
  if (is.null(model_dir)) die("--model-dir is required with --backend=hf")
  if (!dir.exists(model_dir)) die("model dir not found: ", model_dir)
}

ontology <- NULL
if (!is_unset(opts$ontology)) {
  if (!file.exists(opts$ontology)) die("ontology file not found: ", opts$ontology)
  ontology <- paste(readLines(opts$ontology, warn = FALSE), collapse = " ")
}

show <- function(lbl, value) {
  pretty <- if (is.null(value)) "<none>" else paste(value, collapse = ",")
  if (nchar(pretty) > 70) pretty <- paste0(substr(pretty, 1, 67), "...")
  message(sprintf("  %-18s : %s", lbl, pretty))
}
message("== Gemma phenotyper config ==")
show("spe_file",    opts$spe_file)
show("out_dir",     opts$out_dir)
show("backend",     backend)
show("prompt_style", prompt_style)
show("mode",        mode)
show("cluster_col", cluster_col)
show("markers",     if (is.null(markers)) "<all>" else markers)
show("tissue",      tissue)
show("label",       label)
if (backend == "ollama") {
  show("ollama_host",  ollama_host)
  show("ollama_model", ollama_model)
  show("temperature",  temperature)
  show("timeout",      timeout)
} else {
  show("model_dir",    model_dir)
  show("base_model",   opts$base_model)
  show("load_in_8bit", load_in_8bit)
  show("python",       python_bin)
}
# Backend-independent, so echoed for both.
show("harmonize",   harmonize)
show("retry",       retry)
if (prompt_style == "zero-shot") { show("min_z", min_z); show("marker_glossary", marker_glossary) }
show("vote_samples", vote_samples)
if (vote_samples > 1) {
  show("vote_temperature",   vote_temperature)
  show("vote_min_agreement", vote_min_agreement)
  show("vote_override_floor", vote_override_floor)
}
show("export_loom", export_loom)
message("=============================")

## ---------------------------------------------------------------------------
## Load pipeline dependencies (deferred so --help / arg errors stay clean)
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  suppressWarnings({
    library(SpatialExperiment)
    library(SingleCellExperiment)
    library(jsonlite)
    library(curl)
  })
})

## ---------------------------------------------------------------------------
## Source pipeline code
## ---------------------------------------------------------------------------

for (f in c("gemma-utils.R", "RunGemmaPhenotyper.R", "vectra-export.R", "loom-export.R")) {
  path <- file.path(script_dir, f)
  if (!file.exists(path))
    die(f, " not found under '", script_dir, "' (set GEMMA_DIR to override)")
  source(path)
}

infer_script <- file.path(script_dir, "infer.py")
if (backend == "local" && !file.exists(infer_script))
  die("infer.py not found under '", script_dir, "' (set GEMMA_DIR to override)")

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

RunGemmaPhenotyper(
  spe_file     = opts$spe_file,
  out_dir      = opts$out_dir,
  backend      = backend,
  prompt_style = prompt_style,
  label        = label,
  cluster_col  = cluster_col,
  marker_cols  = markers,
  mode         = mode,
  ollama_host  = ollama_host,
  ollama_model = ollama_model,
  temperature  = temperature,
  timeout      = timeout,
  retry        = retry,
  min_z        = min_z,
  marker_glossary = marker_glossary,
  vote_samples = vote_samples,
  vote_temperature = vote_temperature,
  vote_min_agreement = vote_min_agreement,
  vote_override_floor = vote_override_floor,
  harmonize    = harmonize,
  tissue       = tissue,
  model_dir    = model_dir,
  base_model   = opts$base_model,
  load_in_8bit = load_in_8bit,
  ontology     = ontology,
  python_bin   = python_bin,
  infer_script = infer_script,
  export_loom  = export_loom
)
