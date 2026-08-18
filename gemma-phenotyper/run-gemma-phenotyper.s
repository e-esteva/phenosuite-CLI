#!/bin/bash
#
# run-gemma-phenotyper.s — SLURM array task for the Gemma mIF phenotyper.
#
# Submitted by gemma-phenotyper-meta.s with --array=1-${batch_size}. Each task
# extracts its own row (SLURM_ARRAY_TASK_ID) from every batch-input file and
# calls run-gemma-phenotyper.R with named flags.
#
#SBATCH --error=GemmaPhenotyper_%A_%a.err
#SBATCH --out=GemmaPhenotyper_%A_%a.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GEMMA_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/gemma-phenotyper-config.txt"

module load r/4.1.2

# ---------------------------------------------------------------------------
# Extract row N from a line-aligned batch-input file. Prints an empty string
# (which run-gemma-phenotyper.R treats as 'not set') if the file has fewer than
# N lines, and logs a warning so the cause is visible in the SLURM log.
# ---------------------------------------------------------------------------
extract_row() {
    local file="$1"
    local row="$2"
    local label="$3"

    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "run-gemma-phenotyper.s: ${label} list not found: ${file}" >&2
        echo ""
        return
    fi
    local total
    total=$(wc -l < "${file}")
    if (( row > total )); then
        echo "run-gemma-phenotyper.s: ${label} has only ${total} rows, task ${row} will use empty value" >&2
        echo ""
        return
    fi
    sed -n "${row}p" "${file}"
}

TASK="${SLURM_ARRAY_TASK_ID}"

spe_file_tmp=$(extract_row    "${spe_file}"     "${TASK}" "spe_files")
out_dir_tmp=$(extract_row     "${out_dir}"      "${TASK}" "out_dirs")
label_tmp=$(extract_row       "${labels}"       "${TASK}" "labels")
cluster_col_tmp=$(extract_row "${cluster_cols}" "${TASK}" "cluster_cols")
tissue_tmp=$(extract_row      "${tissues}"      "${TASK}" "tissues")

if [[ -z "${spe_file_tmp}" ]]; then
    echo "run-gemma-phenotyper.s: no spe file for task ${TASK}; aborting" >&2
    exit 1
fi
if [[ -z "${out_dir_tmp}" ]]; then
    echo "run-gemma-phenotyper.s: no out_dir for task ${TASK}; aborting" >&2
    exit 1
fi

mkdir -p "${out_dir_tmp}"

echo "== run-gemma-phenotyper.s task ${TASK} =="
echo "  spe_file     : ${spe_file_tmp}"
echo "  out_dir      : ${out_dir_tmp}"
echo "  label        : ${label_tmp:-<basename>}"
echo "  cluster_col  : ${cluster_col_tmp:-<none>}"
echo "  tissue       : ${tissue_tmp:-<none>}"
echo "  backend      : ${backend}"
echo "  prompt_style : ${prompt_style}"
echo "  mode         : ${mode}"
echo "  harmonize    : ${harmonize}"
echo "  vote_samples : ${vote_samples:-0}"
if [[ "${backend}" == "ollama" ]]; then
echo "  ollama_host  : ${ollama_host}"
echo "  ollama_model : ${ollama_model}"
echo "  temperature  : ${temperature}"
else
echo "  model_dir    : ${model_dir}"
echo "  base_model   : ${base_model}"
echo "  load_in_8bit : ${load_in_8bit}"
echo "  gemma_python : ${gemma_python}"
echo "  hf_home      : ${hf_home:-<default: \$HOME/.cache/huggingface>}"
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null \
    || echo "  (no GPU visible — a large checkpoint will fall back to CPU and be very slow)"
fi

# Backend-specific flags: only pass the set the chosen backend actually reads,
# so an unset value for the other backend can never be misread as a real one.
if [[ "${backend}" == "ollama" ]]; then
    backend_flags=(
        --ollama-host="${ollama_host}"
        --ollama-model="${ollama_model}"
        --temperature="${temperature:-0.2}"
    )
else
    export GEMMA_PYTHON="${gemma_python:-python3}"
    # Point HF at project storage before anything imports transformers. Without
    # this the cache defaults to $HOME/.cache/huggingface, which on most
    # clusters is quota'd well below a 12B checkpoint.
    if [[ -n "${hf_home:-}" ]]; then
        export HF_HOME="${hf_home}"
    fi
    # Keep HF fully offline on compute nodes: without these, transformers will
    # try to reach huggingface.co to revalidate the cache and hang or fail on a
    # node with no outbound network. Weights must already be on disk.
    export HF_HUB_OFFLINE="${HF_HUB_OFFLINE:-1}"
    export TRANSFORMERS_OFFLINE="${TRANSFORMERS_OFFLINE:-1}"
    backend_flags=(
        --model-dir="${model_dir}"
        --base-model="${base_model}"
        --load-in-8bit="${load_in_8bit:-false}"
        --ontology="${ontology:-}"
        --python="${gemma_python:-python3}"
    )
fi

Rscript "${SCRIPT_DIR}/run-gemma-phenotyper.R" \
    --spe-file="${spe_file_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --label="${label_tmp}" \
    --cluster-col="${cluster_col_tmp}" \
    --tissue="${tissue_tmp}" \
    --backend="${backend}" \
    --prompt-style="${prompt_style:-zero-shot}" \
    --mode="${mode:-cluster}" \
    --harmonize="${harmonize:-true}" \
    --min-z="${min_z:-0.5}" \
    --marker-glossary="${marker_glossary:-}" \
    --vote-samples="${vote_samples:-0}" \
    --vote-temperature="${vote_temperature:-0.7}" \
    --vote-min-agreement="${vote_min_agreement:-0.6}" \
    --vote-override-floor="${vote_override_floor:-false}" \
    --markers="${markers:-}" \
    --export-loom="${export_loom:-false}" \
    "${backend_flags[@]}"
