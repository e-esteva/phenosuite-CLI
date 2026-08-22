#!/bin/bash
#
# run-segmentation.s — SLURM array task for the segmentation pipeline.
#
# Submitted by segmentation-meta.s with --array=1-${batch_size}. Each task
# extracts its own row (SLURM_ARRAY_TASK_ID) from every batch-input file and
# calls run_segmentation.py with named flags.
#
#SBATCH --error=segmentation_%A_%a.err
#SBATCH --out=segmentation_%A_%a.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SEGMENTATION_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/segmentation-config.txt"

# conda's activation machinery is not written to survive `set -u`: its hook
# and activate.d scripts reference variables an interactive shell would have
# but a batch job never sets, so activation aborts the task outright with
#   slurm_script: line NNN: PS1: unbound variable
# `source` is also a POSIX special builtin, so under `set -e` a missing
# `activate` kills the shell before any `|| fallback` can run. Neither is our
# code, so the mechanism is checked first and nounset is lifted just around
# the call. Errexit stays on.
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
        echo "run-segmentation.s: conda is not on PATH; using the inherited environment" >&2
        ok=1
    fi
    set -u
    if (( ok == 0 )); then
        echo "run-segmentation.s: could not activate conda env '${env_name}'; using the inherited environment" >&2
    fi
}

activate_conda segmentation

if ! command -v python >/dev/null 2>&1; then
    echo "run-segmentation.s: python not found after activation." >&2
    echo "  Create the env (conda env create -f environment.yml) or activate it before submitting" >&2
    echo "  — segmentation-meta.s submits with --export=ALL, so the submitting shell carries over." >&2
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
        echo "run-segmentation.s: ${label} list not found: ${file}" >&2
        echo ""; return
    fi
    local total
    total=$(wc -l < "${file}")
    if (( row > total )); then
        echo "run-segmentation.s: ${label} has only ${total} rows; task ${row} uses empty value" >&2
        echo ""; return
    fi
    sed -n "${row}p" "${file}"
}

TASK="${SLURM_ARRAY_TASK_ID}"

image_tmp=$(   extract_row "${image}"            "${TASK}" "image")
out_dir_tmp=$( extract_row "${out_dirs}"         "${TASK}" "out_dirs")
label_tmp=$(   extract_row "${labels}"           "${TASK}" "labels")
method_tmp=$(  extract_row "${methods}"          "${TASK}" "methods")
nuc_tmp=$(     extract_row "${nuclear_channels}" "${TASK}" "nuclear_channels")
mem_tmp=$(     extract_row "${membrane_channels}" "${TASK}" "membrane_channels")

if [[ -z "${image_tmp}" ]]; then
    echo "run-segmentation.s: no image for task ${TASK}; aborting" >&2
    exit 1
fi
if [[ -z "${out_dir_tmp}" ]]; then
    echo "run-segmentation.s: no out_dir for task ${TASK}; aborting" >&2
    exit 1
fi

# Per-image method falls back to the global default_method when NULL/empty.
if [[ -z "${method_tmp}" || "${method_tmp}" == "NULL" ]]; then
    method_tmp="${default_method}"
fi

mkdir -p "${out_dir_tmp}"

echo "== run-segmentation.s task ${TASK} =="
echo "  image             : ${image_tmp}"
echo "  out_dir           : ${out_dir_tmp}"
echo "  label             : ${label_tmp:-<image basename>}"
echo "  method            : ${method_tmp}"
echo "  nuclear_channel   : ${nuc_tmp:-<auto>}"
echo "  membrane_channel  : ${mem_tmp:-<none>}"
echo "  gpu               : ${gpu}"

gpu_upper=$(echo "${gpu}" | tr '[:lower:]' '[:upper:]')
gpu_flag=()
case "${gpu_upper}" in
    TRUE)  gpu_flag=(--gpu) ;;
    FALSE) gpu_flag=(--no-gpu) ;;
    *)     gpu_flag=() ;;  # AUTO -> let run_segmentation.py auto-detect
esac

overlay_flag=()
[[ "$(echo "${export_overlay}" | tr '[:lower:]' '[:upper:]')" == "FALSE" ]] && overlay_flag=(--no-overlay)

centroids_flag=()
[[ "$(echo "${export_centroids}" | tr '[:lower:]' '[:upper:]')" == "FALSE" ]] && centroids_flag=(--no-centroids)

python "${SCRIPT_DIR}/run_segmentation.py" \
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
    --sam-max-side="${sam_max_side}" \
    --resolution="${resolution}" \
    --min-size="${min_size}" \
    --tile-size="${tile_size}" \
    --tile-overlap="${tile_overlap}" \
    --seed="${seed}" \
    "${gpu_flag[@]}" "${overlay_flag[@]}" "${centroids_flag[@]}"

bash "${SCRIPT_DIR}/makeRunLog-batch.sh"
