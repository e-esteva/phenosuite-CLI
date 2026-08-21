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

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PCF_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/pcf-config.txt"

module load r/4.4.0

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

mkdir -p "${out_dir_tmp}"

echo "== run-pcf.s task ${TASK} =="
echo "  vectra_files    : ${vectra_files_tmp}"
echo "  out_dir         : ${out_dir_tmp}"
echo "  label           : ${label_tmp:-<date>}"
echo "  celltypes       : ${celltypes_tmp:-<shared across inputs>}"
echo "  ref_celltype    : ${ref_celltype_tmp:-<first celltype>}"
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
