#!/bin/bash
#
# run-merfish.s — SLURM array task for the MERFISH pipeline.
#
# Submitted by merfish-meta.s with --array=1-${batch_size}. Each task extracts
# its own row (SLURM_ARRAY_TASK_ID) from every batch-input list and calls
# run-merfish.R with named flags built from merfish-config.txt.
#
#SBATCH --error=RunMerfish_%A_%a.err
#SBATCH --out=RunMerfish_%A_%a.out

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MERFISH_DIR="${SCRIPT_DIR}"

source "${SCRIPT_DIR}/merfish-config.txt"

module load r/4.1.2 2>/dev/null || true

# ---------------------------------------------------------------------------
# Extract row N from a line-aligned batch-input file. Prints an empty string
# (treated as 'not set' downstream) when the file is missing or too short, and
# logs a warning so the cause is visible in the SLURM log.
# ---------------------------------------------------------------------------
extract_row() {
    local file="$1" row="$2" label="$3"
    if [[ -z "${file}" || ! -f "${file}" ]]; then
        echo "run-merfish.s: ${label} list not found: ${file}" >&2
        echo ""; return
    fi
    local total
    total=$(wc -l < "${file}")
    if (( row > total )); then
        echo "run-merfish.s: ${label} has only ${total} rows; task ${row} uses empty value" >&2
        echo ""; return
    fi
    sed -n "${row}p" "${file}"
}

# 'NULL' (the manifest pad token) -> empty, so run-merfish.R applies its default.
denull() { [[ "$1" == "NULL" ]] && echo "" || echo "$1"; }

TASK="${SLURM_ARRAY_TASK_ID:?run-merfish.s must run inside a SLURM array (SLURM_ARRAY_TASK_ID unset)}"

expression_file_tmp=$(extract_row "${expression_file}" "${TASK}" "expression_files")
metadata_file_tmp=$(extract_row   "${metadata_file}"   "${TASK}" "metadata_files")
out_dir_tmp=$(extract_row         "${out_dir}"         "${TASK}" "out_dirs")
sample_id_tmp=$(denull "$(extract_row "${sample_id}"   "${TASK}" "sample_ids")")

if [[ -z "${expression_file_tmp}" ]]; then
    echo "run-merfish.s: no expression file for task ${TASK}; aborting" >&2; exit 1
fi
if [[ -z "${metadata_file_tmp}" ]]; then
    echo "run-merfish.s: no metadata file for task ${TASK}; aborting" >&2; exit 1
fi
if [[ -z "${out_dir_tmp}" ]]; then
    echo "run-merfish.s: no out_dir for task ${TASK}; aborting" >&2; exit 1
fi
mkdir -p "${out_dir_tmp}"

# Optional genes-x-cells switch.
transpose_flag=()
[[ "${transpose:-FALSE}" == "TRUE" ]] && transpose_flag=(--transpose)

echo "== run-merfish.s task ${TASK} =="
echo "  sample_id        : ${sample_id_tmp:-<out_dir basename>}"
echo "  expression_file  : ${expression_file_tmp}"
echo "  metadata_file    : ${metadata_file_tmp}"
echo "  out_dir          : ${out_dir_tmp}"
echo "  x_col / y_col    : ${x_col} / ${y_col}"
echo "  norm_method      : ${norm_method}"
echo "  cluster          : k=${cluster_k} res=${cluster_res} ${leiden_objective}"

Rscript "${SCRIPT_DIR}/run-merfish.R" \
    --expression-file="${expression_file_tmp}" \
    --metadata-file="${metadata_file_tmp}" \
    --out-dir="${out_dir_tmp}" \
    ${sample_id_tmp:+--sample-id="${sample_id_tmp}"} \
    "${transpose_flag[@]}" \
    --x-col="${x_col}" \
    --y-col="${y_col}" \
    --area-col="$(denull "${area_col}")" \
    --negctrl-col="$(denull "${negctrl_col}")" \
    --vol-col="$(denull "${vol_col}")" \
    --qc-min-counts="${qc_min_counts}" \
    --qc-max-counts="${qc_max_counts}" \
    --qc-min-genes="${qc_min_genes}" \
    --qc-max-genes="${qc_max_genes}" \
    --qc-min-area="${qc_min_area}" \
    --qc-max-area="${qc_max_area}" \
    --qc-min-density="${qc_min_density}" \
    --qc-max-density="${qc_max_density}" \
    --qc-max-negctrl-ratio="${qc_max_negctrl_ratio}" \
    --norm-method="${norm_method}" \
    --scale-method="${scale_method}" \
    --hvg-method="${hvg_method}" \
    --n-hvg="${n_hvg}" \
    --n-pcs="${n_pcs}" \
    --umap-neighbors="${umap_neighbors}" \
    --umap-min-dist="${umap_min_dist}" \
    --cluster-k="${cluster_k}" \
    --cluster-res="${cluster_res}" \
    --leiden-objective="${leiden_objective}" \
    --nhood-k="${nhood_k}" \
    --svg-k="${svg_k}" \
    --svg-n-top="${svg_n_top}" \
    --export-spe="${export_spe}" \
    --export-figures="${export_figures}" \
    --fig-width="${fig_width}" \
    --fig-height="${fig_height}" \
    --fig-dpi="${fig_dpi}" \
    --seed="${seed}"
