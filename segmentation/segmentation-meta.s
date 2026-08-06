#!/bin/bash
#
# segmentation-meta.s — entry point for the segmentation batch pipeline.
#
# Reads segmentation-config.txt, validates that every batch-input file is
# aligned (has at least batch_size rows), then submits run-segmentation.s
# as a SLURM array job sized 1..batch_size.
#
# Usage:
#   sbatch segmentation-meta.s
#
#SBATCH --time=0-1
#SBATCH --mem=4GB
#SBATCH --error=segmentation-meta_%j.err
#SBATCH --out=segmentation-meta_%j.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export SEGMENTATION_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/segmentation-config.txt"

# ---------------------------------------------------------------------------
# Pre-submit validation: every per-image list must have >= batch_size rows so
# task N never runs on a phantom row.
# ---------------------------------------------------------------------------
check_rows() {
    local file="$1"
    local label="$2"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "segmentation-meta.s: ${label} file not found: ${file}" >&2
        exit 1
    fi
    local total
    total=$(wc -l < "${file}")
    if (( total < batch_size )); then
        echo "segmentation-meta.s: ${label} has ${total} rows but batch_size=${batch_size}." >&2
        echo "  Pad the file to ${batch_size} lines (use 'NULL' for empty rows)," >&2
        echo "  or reduce the primary input list so batch_size matches." >&2
        exit 1
    fi
}

check_rows "${image}"            "image"
check_rows "${out_dirs}"         "out_dirs"
check_rows "${labels}"           "labels"
check_rows "${methods}"          "methods"
check_rows "${nuclear_channels}" "nuclear_channels"
check_rows "${membrane_channels}" "membrane_channels"

echo "segmentation-meta.s: batch_size=${batch_size}, all per-image inputs aligned."
echo "Initiating segmentation at $(date)"

gres_flag=()
if [[ -n "${module1_gres:-}" ]]; then
    gres_flag=(--gres="${module1_gres}")
fi

job_id=$(sbatch \
    --export=ALL,configfile="${configFile}" \
    --mem="${module1_mem}" \
    --time="${module1_time}" \
    --array="1-${batch_size}" \
    --partition="${module1_partition}" \
    "${gres_flag[@]}" \
    "${SCRIPT_DIR}/${module1_Path}" \
    | awk '{print $NF}')

echo "Submitted array job: ${job_id}"
echo "Finished submission at $(date)"
