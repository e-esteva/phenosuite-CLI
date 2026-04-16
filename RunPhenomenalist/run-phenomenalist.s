#!/bin/bash
#
# run-phenomenalist.s — SLURM array task for RunPhenomenalist.
#
# Submitted by phenomenalist-meta.s with --array=1-${batch_size}. Each task
# extracts its own row (SLURM_ARRAY_TASK_ID) from every batch-input file and
# calls run-phenomenalist.R with named flags.
#
#SBATCH --error=RunPhenomenalist_%A_%a.err
#SBATCH --out=RunPhenomenalist_%A_%a.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PHENOMENALIST_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/phenomenalist-config.txt"

module load r/4.1.2

# ---------------------------------------------------------------------------
# Extract row N from a line-aligned batch-input file. Prints an empty string
# (which run-phenomenalist.R treats as 'not set') if the file has fewer than N
# lines, and logs a warning so the cause is visible in the SLURM log.
# ---------------------------------------------------------------------------
extract_row() {
    local file="$1"
    local row="$2"
    local label="$3"

    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "run-phenomenalist.s: ${label} list not found: ${file}" >&2
        echo ""
        return
    fi
    local total
    total=$(wc -l < "${file}")
    if (( row > total )); then
        echo "run-phenomenalist.s: ${label} has only ${total} rows, task ${row} will use empty value" >&2
        echo ""
        return
    fi
    sed -n "${row}p" "${file}"
}

TASK="${SLURM_ARRAY_TASK_ID}"

segmentation_file_tmp=$(extract_row "${segmentation_file}" "${TASK}" "segmentation_files")
failed_markers_tmp=$(extract_row   "${failed_markers}"    "${TASK}" "failed-markers")
nuclear_markers_tmp=$(extract_row  "${nuclear_markers}"   "${TASK}" "nuclear-markers")
out_dir_tmp=$(extract_row          "${out_dir}"           "${TASK}" "out_dirs")
classifier_label_tmp=$(extract_row "${classifier_labels}" "${TASK}" "labels")

if [[ -z "${segmentation_file_tmp}" ]]; then
    echo "run-phenomenalist.s: no segmentation file for task ${TASK}; aborting" >&2
    exit 1
fi
if [[ -z "${out_dir_tmp}" ]]; then
    echo "run-phenomenalist.s: no out_dir for task ${TASK}; aborting" >&2
    exit 1
fi

mkdir -p "${out_dir_tmp}"

echo "== run-phenomenalist.s task ${TASK} =="
echo "  segmentation_file : ${segmentation_file_tmp}"
echo "  failed_markers    : ${failed_markers_tmp:-<none>}"
echo "  nuclear_markers   : ${nuclear_markers_tmp:-<none>}"
echo "  out_dir           : ${out_dir_tmp}"
echo "  classifier_label  : ${classifier_label_tmp:-<none>}"
echo "  clustering_res    : ${clustering_res}"
echo "  max_cells         : ${max_cells}"
echo "  phenotyping_tmpl  : ${phenotyping_template:-<none>}"
echo "  skip_cols         : ${skip_cols:-<auto>}"

Rscript "${SCRIPT_DIR}/run-phenomenalist.R" \
    --segmentation-file="${segmentation_file_tmp}" \
    --failed-markers="${failed_markers_tmp}" \
    --nuclear-markers="${nuclear_markers_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --clustering-res="${clustering_res}" \
    --classifier-label="${classifier_label_tmp}" \
    --max-cells="${max_cells}" \
    --phenotyping-template="${phenotyping_template:-}" \
    --skip-cols="${skip_cols:-}"
