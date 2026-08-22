#!/bin/bash
#
# run-pcf.s — SLURM array task for the PCF pipeline.
#
# Submitted by pcf-meta.s with --array=1-${batch_size}. Each task extracts its
# own row (SLURM_ARRAY_TASK_ID) from every batch-input file and calls run-pcf.R
# with named flags.
#
#SBATCH --error=PCF_%A_%a.err
#SBATCH --out=PCF_%A_%a.out

set -euo pipefail

# SLURM copies the batch script into its spool directory and runs it from
# there, so ${BASH_SOURCE[0]} is /…/spoold/…/slurm_script — not this module's
# directory, and the config file is not beside it. Resolve the real location
# by looking for pcf-config.txt in each candidate that could mean something,
# most authoritative first. PCF_DIR comes first so pcf-meta.s can hand the
# array tasks the answer it already resolved (it submits with --export=ALL).
find_pcf_dir() {
    local c
    for c in "${PCF_DIR:-}" \
             "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" \
             "${SLURM_SUBMIT_DIR:-}" \
             "${SLURM_SUBMIT_DIR:+${SLURM_SUBMIT_DIR}/pcf}" \
             "${PWD}" \
             "${PWD}/pcf"; do
        if [[ -n "${c}" && -f "${c}/pcf-config.txt" ]]; then
            printf '%s\n' "${c}"
            return 0
        fi
    done
    return 1
}

if ! SCRIPT_DIR="$(find_pcf_dir)"; then
    echo "run-pcf.s: cannot locate pcf-config.txt." >&2
    echo "  Searched: \$PCF_DIR, this script's directory, \$SLURM_SUBMIT_DIR (+ its pcf/), \$PWD (+ its pcf/)." >&2
    echo "  SLURM runs a *copy* of the batch script from its spool directory, so the config is" >&2
    echo "  not next to it. Submit from the module directory (cd pcf && sbatch pcf-meta.s)," >&2
    echo "  or export PCF_DIR=/path/to/pcf before submitting." >&2
    exit 1
fi
export PCF_DIR="${SCRIPT_DIR}"

# pcf-config.txt declares its batch-input paths relative to the module
# directory (batch-inputs/…), and makeRunLog-batch.sh sources the config the
# same way. A SLURM task starts in the submit directory, which is not
# necessarily that one, so anchor here before anything reads those paths.
cd "${SCRIPT_DIR}"

source "${SCRIPT_DIR}/pcf-config.txt"

# ---------------------------------------------------------------------------
# R provisioning. Both keys are optional and both may be left empty.
#
# This module ships its own conda env (environment.yml), and a `module load`
# would put the site R ahead of that env's R on PATH — hiding the exact
# spatstat / ggpubr install the env exists to provide. So nothing is loaded
# unless it is asked for:
#   conda_env  — env to activate in each array task (default: runpcf)
#   r_module   — site R module for clusters that keep R in modules instead
#
# Leaving both empty inherits whatever the submitting shell had: pcf-meta.s
# submits with --export=ALL, so an env activated before `sbatch` carries into
# every task on its own.
# ---------------------------------------------------------------------------
# `source` is a POSIX special builtin: under `set -e` a failed `source activate`
# aborts the task immediately, before any `|| fallback` can run. So the
# mechanism is checked before it is used, and a missing conda degrades to a
# warning instead of a silent exit 1.
# conda's activation machinery is not written to survive `set -u`: its hook
# and activate.d scripts reference variables an interactive shell would have
# but a batch job never sets, so activation aborts the task outright with
#   slurm_script: line NNN: PS1: unbound variable
# That code is not ours to fix, so nounset is lifted just around the call and
# PS1 is given an empty default. Errexit stays on.
activate_conda() {
    local env_name="$1"
    : "${PS1:=}"
    set +u
    local ok=0
    if command -v conda >/dev/null 2>&1; then
        eval "$(conda shell.bash hook)" && conda activate "${env_name}" && ok=1
    elif command -v activate >/dev/null 2>&1; then
        source activate "${env_name}" && ok=1
    else
        echo "run-pcf.s: conda_env='${env_name}' is set but conda is not on PATH; using the inherited environment" >&2
        ok=1
    fi
    set -u
    if (( ok == 0 )); then
        echo "run-pcf.s: could not activate conda env '${env_name}'; using the inherited environment" >&2
    fi
}

if [[ -n "${conda_env:-}" && "${conda_env}" != "NULL" ]]; then
    activate_conda "${conda_env}"
fi

# `module` is usually a shell function from /etc/profile.d, which a
# non-interactive job shell may not have sourced. An r_module that was asked
# for but cannot be loaded is an error, not something to run past quietly.
if [[ -n "${r_module:-}" && "${r_module}" != "NULL" ]]; then
    if [[ "$(type -t module || true)" == "" ]]; then
        echo "run-pcf.s: r_module='${r_module}' is set but no 'module' command is available here." >&2
        echo "  Clear r_module= in pcf-config.txt and use conda_env= instead, or source the" >&2
        echo "  module system (e.g. /etc/profile.d/modules.sh) before submitting." >&2
        exit 1
    fi
    # Module implementations are shell functions with the same nounset
    # problem as conda's hook; lift -u around the load for the same reason.
    set +u
    module load "${r_module}"
    set -u
fi

if ! command -v Rscript >/dev/null 2>&1; then
    echo "run-pcf.s: Rscript not found." >&2
    echo "  Set conda_env= (e.g. runpcf) or r_module= (e.g. r/4.4.0) in pcf-config.txt," >&2
    echo "  or activate the environment before submitting — pcf-meta.s submits with" >&2
    echo "  --export=ALL, so the submitting shell's environment carries into each task." >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Extract row N from a line-aligned batch-input file.
# Prints empty string (treated as 'not set') when the file is shorter than N.
# ---------------------------------------------------------------------------
extract_row() {
    local file="$1"
    local row="$2"
    local label="$3"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "run-pcf.s: ${label} list not found: ${file}" >&2
        echo ""; return
    fi
    local total
    total=$(wc -l < "${file}")
    if (( row > total )); then
        echo "run-pcf.s: ${label} has only ${total} rows; task ${row} uses empty value" >&2
        echo ""; return
    fi
    sed -n "${row}p" "${file}"
}

TASK="${SLURM_ARRAY_TASK_ID}"

vectra_files_tmp=$( extract_row "${vectra_files}"  "${TASK}" "vectra_files")
out_dir_tmp=$(      extract_row "${out_dirs}"      "${TASK}" "out_dirs")
label_tmp=$(        extract_row "${labels}"        "${TASK}" "labels")
celltypes_tmp=$(    extract_row "${celltypes}"     "${TASK}" "celltypes")
ref_celltype_tmp=$( extract_row "${ref_celltypes}" "${TASK}" "ref_celltypes")

if [[ -z "${vectra_files_tmp}" ]]; then
    echo "run-pcf.s: no vectra_files for task ${TASK}; aborting" >&2
    exit 1
fi
if [[ -z "${out_dir_tmp}" ]]; then
    echo "run-pcf.s: no out_dir for task ${TASK}; aborting" >&2
    exit 1
fi

# Values are echoed with [brackets] *before* anything acts on them: a stray
# CR, a trailing space, or a row picked from a misaligned file is invisible in
# a bare path, and mkdir failing under `set -e` would otherwise kill the task
# before the diagnostics below ever printed.
echo "== run-pcf.s task ${TASK} =="
echo "  config dir      : [${SCRIPT_DIR}]"
echo "  vectra_files    : [${vectra_files_tmp}]"
echo "  out_dir         : [${out_dir_tmp}]"
echo "  label           : [${label_tmp}]"
echo "  celltypes       : [${celltypes_tmp}]"
echo "  ref_celltype    : [${ref_celltype_tmp}]"

if ! mkdir -p "${out_dir_tmp}"; then
    echo "run-pcf.s: could not create out_dir [${out_dir_tmp}]" >&2
    echo "  Read from row ${TASK} of ${out_dirs}. If that is not the path you expect," >&2
    echo "  the batch-input files are misaligned — check 'wc -l batch-inputs/*.txt'" >&2
    echo "  and 'sed -n \"${TASK}p\" ${out_dirs} | cat -A' for stray CR/whitespace." >&2
    exit 1
fi

echo "== resolved parameters =="
echo "  radius          : ${radius}"
echo "  resolution      : ${resolution}"
echo "  count_threshold : ${count_threshold}"
echo "  min_count       : ${min_count}"
echo "  phenotype_col   : ${phenotype_col}"
echo "  make_plots      : ${make_plots}"

Rscript "${SCRIPT_DIR}/run-pcf.R" \
    --vectra-files="${vectra_files_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --label="${label_tmp}" \
    --celltypes="${celltypes_tmp}" \
    --ref-celltype="${ref_celltype_tmp}" \
    --radius="${radius}" \
    --resolution="${resolution}" \
    --count-threshold="${count_threshold}" \
    --min-count="${min_count}" \
    --phenotype-col="${phenotype_col}" \
    $( [[ "${make_plots}" == "FALSE" ]] && echo "--no-plots" )

bash "${SCRIPT_DIR}/makeRunLog-batch.sh"
