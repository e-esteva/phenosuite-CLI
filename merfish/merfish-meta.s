#!/bin/bash
#
# merfish-meta.s — entry point for the MERFISH batch pipeline.
#
# Reads merfish-config.txt, validates that the batch-input lists line up, and
# submits run-merfish.s as a SLURM array job sized 1..batch_size.
#
# Usage:
#   sbatch merfish-meta.s
#
#SBATCH --time=0-1
#SBATCH --mem=4GB
#SBATCH --error=RunMerfish-meta_%j.err
#SBATCH --out=RunMerfish-meta_%j.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MERFISH_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/merfish-config.txt"

# ---------------------------------------------------------------------------
# Pre-submit validation: every per-sample list must have >= batch_size rows so
# task N never runs on a phantom row.
# ---------------------------------------------------------------------------
check_rows() {
    local file="$1" label="$2"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "merfish-meta.s: ${label} file not found: ${file}" >&2; exit 1
    fi
    local total
    total=$(wc -l < "${file}")
    if (( total < batch_size )); then
        echo "merfish-meta.s: ${label} has ${total} rows but batch_size=${batch_size}." >&2
        echo "  Pad it to ${batch_size} lines (use 'NULL' for empty rows), or shrink" >&2
        echo "  the primary expression list so batch_size matches." >&2
        exit 1
    fi
}

check_rows "${expression_file}" "expression_file"
check_rows "${metadata_file}"   "metadata_file"
check_rows "${out_dir}"         "out_dir"
# sample_id is optional; validate only if a real path is configured.
if [[ -n "${sample_id:-}" && -f "${sample_id}" ]]; then
    check_rows "${sample_id}" "sample_id"
fi

echo "merfish-meta.s: batch_size=${batch_size}, all per-sample inputs aligned."
echo "Initiating MERFISH batch at $(date)"

submit_merfish() {
    sbatch \
        --export=ALL,configfile="${configFile}" \
        --cpus-per-task="${module1_cpus}" \
        --mem="${module1_mem}" \
        --time="${module1_time}" \
        --array="1-${batch_size}" \
        --partition="${module1_partition}" \
        "${SCRIPT_DIR}/${module1_Path}" \
        | awk '{print $NF}'
}

mod1_job=$(submit_merfish)
echo "Submitted array job: ${mod1_job}  (tasks 1-${batch_size})"
echo "Finished submission at $(date)"
