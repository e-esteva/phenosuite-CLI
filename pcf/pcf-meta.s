#!/bin/bash
#
# pcf-meta.s — entry point for the PCF batch pipeline.
#
# Reads pcf-config.txt, validates that every batch-input file is aligned (has
# at least batch_size rows), then submits run-pcf.s as a SLURM array job
# sized 1..batch_size.
#
# Usage:
#   sbatch pcf-meta.s
#
#SBATCH --time=0-1
#SBATCH --mem=4GB
#SBATCH --error=PCF-meta_%j.err
#SBATCH --out=PCF-meta_%j.out

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
    echo "pcf-meta.s: cannot locate pcf-config.txt." >&2
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
# Pre-submit validation
# ---------------------------------------------------------------------------
check_rows() {
    local file="$1"
    local label="$2"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "pcf-meta.s: ${label} file not found: ${file}" >&2
        exit 1
    fi
    local total
    total=$(wc -l < "${file}")
    if (( total < batch_size )); then
        echo "pcf-meta.s: ${label} has ${total} rows but batch_size=${batch_size}." >&2
        echo "  Pad the file to ${batch_size} lines (use 'NULL' for empty rows)," >&2
        echo "  or reduce the primary input list so batch_size matches." >&2
        exit 1
    fi
}

check_rows "${vectra_files}"  "vectra_files"
check_rows "${out_dirs}"      "out_dirs"
check_rows "${labels}"        "labels"
check_rows "${celltypes}"     "celltypes"
check_rows "${ref_celltypes}" "ref_celltypes"

if (( count_threshold < 1 )); then
    echo "pcf-meta.s: count_threshold must be >= 1 (got ${count_threshold})." >&2
    exit 1
fi

echo "pcf-meta.s: batch_size=${batch_size}, all per-run inputs aligned."
echo "Initiating PCF at $(date)"

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
