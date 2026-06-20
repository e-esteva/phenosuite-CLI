#!/bin/bash
#
# neighborhood_analysis-meta.s — entry point for the NeighborhoodR batch pipeline.
#
# Reads neighborhood_analysis-config.txt, validates that every batch-input file
# is aligned (has at least batch_size rows), then submits run-neighborhood_analysis.s
# as a SLURM array job sized 1..batch_size.
#
# Usage:
#   sbatch neighborhood_analysis-meta.s
#
#SBATCH --time=0-1
#SBATCH --mem=4GB
#SBATCH --error=NeighborhoodR-meta_%j.err
#SBATCH --out=NeighborhoodR-meta_%j.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NEIGHBORHOODR_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/neighborhood_analysis-config.txt"

# ---------------------------------------------------------------------------
# Pre-submit validation
# ---------------------------------------------------------------------------
check_rows() {
    local file="$1"
    local label="$2"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "neighborhood_analysis-meta.s: ${label} file not found: ${file}" >&2
        exit 1
    fi
    local total
    total=$(wc -l < "${file}")
    if (( total < batch_size )); then
        echo "neighborhood_analysis-meta.s: ${label} has ${total} rows but batch_size=${batch_size}." >&2
        echo "  Pad the file to ${batch_size} lines (use 'NULL' for empty rows)," >&2
        echo "  or reduce the primary input list so batch_size matches." >&2
        exit 1
    fi
}

check_rows "${rds_files}"      "rds_files"
check_rows "${celltype_cols}"  "celltype_cols"
check_rows "${out_dirs}"       "out_dirs"
check_rows "${labels}"         "labels"
check_rows "${condition_maps}" "condition_maps"

echo "neighborhood_analysis-meta.s: batch_size=${batch_size}, all per-run inputs aligned."
echo "Initiating NeighborhoodR at $(date)"

job_id=$(sbatch \
    --export=ALL,configfile="${configFile}" \
    --mem="${module1_mem}" \
    --time="${module1_time}" \
    --array="1-${batch_size}" \
    --partition="${module1_partition}" \
    "${SCRIPT_DIR}/${module1_Path}" \
    | awk '{print $NF}')

echo "Submitted array job: ${job_id}"
echo "Finished submission at $(date)"
