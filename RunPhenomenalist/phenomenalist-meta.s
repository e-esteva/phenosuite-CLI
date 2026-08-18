#!/bin/bash
#
# phenomenalist-meta.s — entry point for the RunPhenomenalist batch pipeline.
#
# Reads phenomenalist-config.txt, validates the batch-inputs line-up, and
# submits run-phenomenalist.s as a SLURM array job sized 1..batch_size.
#
#SBATCH --time=0-1
#SBATCH --mem=10GB
#SBATCH --error=RunPhenomenalist-meta_%j.err
#SBATCH --out=RunPhenomenalist-meta_%j.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PHENOMENALIST_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/phenomenalist-config.txt"

# ---------------------------------------------------------------------------
# Pre-submit validation: make sure every per-sample batch-input file has at
# least ${batch_size} rows, so task N never runs on a phantom row.
# ---------------------------------------------------------------------------
check_rows() {
    local file="$1"
    local label="$2"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "phenomenalist-meta.s: ${label} file not found: ${file}" >&2
        exit 1
    fi
    local total
    total=$(wc -l < "${file}")
    if (( total < batch_size )); then
        echo "phenomenalist-meta.s: ${label} has ${total} rows but batch_size=${batch_size}." >&2
        echo "  Either pad the file to ${batch_size} lines (use 'NULL' for empty rows)," >&2
        echo "  or reduce the primary input list so batch_size matches." >&2
        exit 1
    fi
}

check_rows "${segmentation_file}"  "segmentation_file"
check_rows "${failed_markers}"     "failed_markers"
check_rows "${nuclear_markers}"    "nuclear_markers"
check_rows "${out_dir}"            "out_dir"
check_rows "${classifier_labels}"  "classifier_labels"

echo "phenomenalist-meta.s: batch_size=${batch_size}, all per-sample inputs aligned."
echo "Initiating phenomenalist at $(date)"

submit_phenomenalist() {
    sbatch \
        --export=ALL,configfile="${configFile}" \
        --mem="${module1_mem}" \
        --time="${module1_time}" \
        --array="1-${batch_size}" \
        --partition="${module1_partition}" \
        "${SCRIPT_DIR}/${module1_Path}" \
        | awk '{print $NF}'
}

mod1_job=$(submit_phenomenalist)
echo "Submitted array job: ${mod1_job}"
echo "Finished submission at $(date)"
