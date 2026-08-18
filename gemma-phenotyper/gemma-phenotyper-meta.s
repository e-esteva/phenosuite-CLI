#!/bin/bash
#
# gemma-phenotyper-meta.s — entry point for the Gemma phenotyper batch pipeline.
#
# Reads gemma-phenotyper-config.txt, validates the batch-inputs line-up and the
# selected backend, and submits run-gemma-phenotyper.s as a SLURM array job
# sized 1..batch_size.
#
#SBATCH --time=0-1
#SBATCH --mem=4GB
#SBATCH --error=GemmaPhenotyper-meta_%j.err
#SBATCH --out=GemmaPhenotyper-meta_%j.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export GEMMA_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/gemma-phenotyper-config.txt"

# ---------------------------------------------------------------------------
# Pre-submit validation: make sure every per-sample batch-input file has at
# least ${batch_size} rows, so task N never runs on a phantom row.
# ---------------------------------------------------------------------------
check_rows() {
    local file="$1"
    local label="$2"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "gemma-phenotyper-meta.s: ${label} file not found: ${file}" >&2
        exit 1
    fi
    local total
    total=$(wc -l < "${file}")
    if (( total < batch_size )); then
        echo "gemma-phenotyper-meta.s: ${label} has ${total} rows but batch_size=${batch_size}." >&2
        echo "  Either pad the file to ${batch_size} lines (use 'NULL' for empty rows)," >&2
        echo "  or reduce the primary input list so batch_size matches." >&2
        exit 1
    fi
}

check_rows "${spe_file}"     "spe_file"
check_rows "${out_dir}"      "out_dir"
check_rows "${labels}"       "labels"
check_rows "${cluster_cols}" "cluster_cols"
check_rows "${tissues}"      "tissues"

# ---------------------------------------------------------------------------
# Backend validation — fail here rather than in ${batch_size} array tasks.
# ---------------------------------------------------------------------------
case "${backend}" in
  hf)
    if [[ -z "${model_dir:-}" ]]; then
        echo "gemma-phenotyper-meta.s: backend=hf but model_dir is unset." >&2
        exit 1
    fi
    if [[ ! -d "${model_dir}" ]]; then
        echo "gemma-phenotyper-meta.s: model_dir does not exist: ${model_dir}" >&2
        exit 1
    fi
    if [[ ! -f "${model_dir}/config.json" && ! -f "${model_dir}/adapter_config.json" ]]; then
        echo "gemma-phenotyper-meta.s: ${model_dir} has neither config.json nor" >&2
        echo "  adapter_config.json — not a recognisable checkpoint or adapter." >&2
        exit 1
    fi
    # Weights must be on disk: compute nodes are typically offline, and the
    # array tasks run with HF_HUB_OFFLINE=1 so nothing will be downloaded.
    # Check each glob separately: `ls a/*.x a/*.y` exits non-zero when EITHER
    # pattern misses, which would reject a normal safetensors-only checkpoint.
    shopt -s nullglob
    weight_files=("${model_dir}"/*.safetensors "${model_dir}"/*.bin)
    shopt -u nullglob
    if (( ${#weight_files[@]} == 0 )); then
        if [[ -f "${model_dir}/adapter_config.json" ]]; then
            echo "gemma-phenotyper-meta.s: NOTE — ${model_dir} is a LoRA adapter; the base" >&2
            echo "  model '${base_model}' must already be in the HF cache on the compute" >&2
            echo "  nodes, since they run offline." >&2
        else
            echo "gemma-phenotyper-meta.s: no *.safetensors / *.bin found in ${model_dir}." >&2
            echo "  Pre-download the weights on a login node before submitting." >&2
            exit 1
        fi
    fi
    if [[ -n "${hf_home:-}" ]]; then
        if [[ ! -d "${hf_home}" ]]; then
            echo "gemma-phenotyper-meta.s: WARNING — hf_home does not exist yet: ${hf_home}" >&2
            echo "  It will be created on first use, but nothing is cached there now." >&2
        fi
        # A LoRA adapter needs its base model resolvable from the cache, and the
        # tasks run offline — so check now rather than failing in every task.
        if [[ -f "${model_dir}/adapter_config.json" ]]; then
            base_slug="models--$(echo "${base_model}" | tr '/' '-' | sed 's/-/--/')"
            if ! ls -d "${hf_home}/hub/models--"* >/dev/null 2>&1; then
                echo "gemma-phenotyper-meta.s: WARNING — model_dir is a LoRA adapter and" >&2
                echo "  ${hf_home}/hub contains no cached models. The base model" >&2
                echo "  '${base_model}' must be downloaded there before submitting," >&2
                echo "  because the array tasks run with HF_HUB_OFFLINE=1." >&2
            fi
        fi
    fi
    if [[ -z "${gemma_python:-}" ]]; then
        echo "gemma-phenotyper-meta.s: backend=hf but gemma_python is unset." >&2
        exit 1
    fi
    if ! "${gemma_python}" -c "import torch, transformers" 2>/dev/null; then
        echo "gemma-phenotyper-meta.s: WARNING — '${gemma_python}' could not import" >&2
        echo "  torch/transformers on the submit host. If the compute nodes use a" >&2
        echo "  different env (module load / conda activate inside the task), ignore" >&2
        echo "  this. Otherwise the array tasks will fail. Submitting anyway." >&2
    fi
    ;;
  ollama)
    if [[ -z "${ollama_model:-}" ]]; then
        echo "gemma-phenotyper-meta.s: backend=ollama but ollama_model is unset." >&2
        exit 1
    fi
    if [[ -z "${ollama_host:-}" ]]; then
        echo "gemma-phenotyper-meta.s: backend=ollama but ollama_host is unset." >&2
        exit 1
    fi
    # Reachability is checked from the submit host. Compute nodes may have
    # different network access, so this is a warning, not a hard failure.
    if command -v curl >/dev/null 2>&1; then
        if curl -sf --max-time 10 "${ollama_host%/}/api/tags" >/dev/null 2>&1; then
            echo "gemma-phenotyper-meta.s: Ollama reachable at ${ollama_host} (from submit host)."
        else
            echo "gemma-phenotyper-meta.s: WARNING — could not reach Ollama at ${ollama_host}" >&2
            echo "  from the submit host. Ensure it is reachable from the COMPUTE nodes," >&2
            echo "  which is what actually matters. Submitting anyway." >&2
        fi
    fi
    ;;
  *)
    echo "gemma-phenotyper-meta.s: backend must be 'hf' or 'ollama' (got '${backend}')" >&2
    exit 1
    ;;
esac

echo "gemma-phenotyper-meta.s: batch_size=${batch_size}, all per-sample inputs aligned."
echo "gemma-phenotyper-meta.s: backend=${backend}, prompt_style=${prompt_style:-zero-shot}"
echo "Initiating gemma-phenotyper at $(date)"

submit_gemma() {
    # --gres is only meaningful for backend=hf; leave module1_gres empty for
    # ollama so the job is not queued behind GPU availability it never uses.
    # Built as a plain string rather than an array: under `set -u`, bash 3.2
    # (still common on login nodes) treats "${empty_array[@]}" as unbound.
    local gres_arg=""
    if [[ -n "${module1_gres:-}" ]]; then
        gres_arg="--gres=${module1_gres}"
    fi
    sbatch \
        --export=ALL,configfile="${configFile}" \
        --mem="${module1_mem}" \
        --time="${module1_time}" \
        --array="1-${batch_size}" \
        --partition="${module1_partition}" \
        ${gres_arg} \
        "${SCRIPT_DIR}/${module1_Path}" \
        | awk '{print $NF}'
}

mod1_job=$(submit_gemma)
echo "Submitted array job: ${mod1_job}"
echo "Finished submission at $(date)"
