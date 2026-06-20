#!/bin/bash
#
# run-neighborhood_analysis.s — SLURM array task for NeighborhoodR.
#
# Submitted by neighborhood_analysis-meta.s with --array=1-${batch_size}.
# Each task extracts its own row (SLURM_ARRAY_TASK_ID) from every batch-input
# file and calls run-neighborhood_analysis.R with named flags.
#
#SBATCH --error=NeighborhoodR_%A_%a.err
#SBATCH --out=NeighborhoodR_%A_%a.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export NEIGHBORHOODR_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/neighborhood_analysis-config.txt"

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
        echo "run-neighborhood_analysis.s: ${label} list not found: ${file}" >&2
        echo ""; return
    fi
    local total
    total=$(wc -l < "${file}")
    if (( row > total )); then
        echo "run-neighborhood_analysis.s: ${label} has only ${total} rows; task ${row} uses empty value" >&2
        echo ""; return
    fi
    sed -n "${row}p" "${file}"
}

TASK="${SLURM_ARRAY_TASK_ID}"

rds_files_tmp=$(    extract_row "${rds_files}"      "${TASK}" "rds_files")
celltype_col_tmp=$( extract_row "${celltype_cols}"  "${TASK}" "celltype_cols")
out_dir_tmp=$(      extract_row "${out_dirs}"        "${TASK}" "out_dirs")
label_tmp=$(        extract_row "${labels}"          "${TASK}" "labels")
cond_map_tmp=$(     extract_row "${condition_maps}"  "${TASK}" "condition_maps")

if [[ -z "${rds_files_tmp}" ]]; then
    echo "run-neighborhood_analysis.s: no rds_files for task ${TASK}; aborting" >&2
    exit 1
fi
if [[ -z "${out_dir_tmp}" ]]; then
    echo "run-neighborhood_analysis.s: no out_dir for task ${TASK}; aborting" >&2
    exit 1
fi

mkdir -p "${out_dir_tmp}"

echo "== run-neighborhood_analysis.s task ${TASK} =="
echo "  rds_files     : ${rds_files_tmp}"
echo "  celltype_col  : ${celltype_col_tmp:-<auto>}"
echo "  out_dir       : ${out_dir_tmp}"
echo "  label         : ${label_tmp:-<date>}"
echo "  k1            : ${k1}"
echo "  k2            : ${k2:-<sweep>}"
echo "  k2_min/max    : ${k2_min} – ${k2_max:-auto}"
echo "  loo_mode      : ${loo_mode}"
echo "  loo_n         : ${loo_n}"
echo "  agg_fn        : ${agg_fn}"
echo "  condition_col : ${condition_col:-<none>}"
echo "  condition_map : ${cond_map_tmp:-<none>}"
echo "  seed          : ${seed}"
echo "  make_plots    : ${make_plots}"

Rscript "${SCRIPT_DIR}/run-neighborhood_analysis.R" \
    --rds-files="${rds_files_tmp}" \
    --celltype-col="${celltype_col_tmp}" \
    --out-dir="${out_dir_tmp}" \
    --label="${label_tmp}" \
    --k1="${k1}" \
    --k2="${k2}" \
    --k2-min="${k2_min}" \
    --k2-max="${k2_max}" \
    --loo-mode="${loo_mode}" \
    --loo-n="${loo_n}" \
    --agg-fn="${agg_fn}" \
    --condition-col="${condition_col}" \
    --condition-map="${cond_map_tmp}" \
    --seed="${seed}" \
    $( [[ "${make_plots}" == "FALSE" ]] && echo "--no-plots" )

bash "${SCRIPT_DIR}/makeRunLog-batch.sh"
