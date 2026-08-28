#!/bin/bash
#SBATCH --error=pwlo_%j.err
#SBATCH --out=pwlo_%j.out

# Activate the shipped conda env. Provision once with:
#   conda env create -f environment.yml
source activate spatial-dynamics

source config-spatial_dynamics.txt

TASK="${SLURM_ARRAY_TASK_ID}"

# Batch-input rows are bare, unquoted paths — do NOT wrap them in quotes in
# the .txt files. Command substitution does not re-parse quotes, so a quoted
# row arrives with the " characters literally inside the filename.
#
# Spaces in a path ("…/CODEX paper revisions/…") are handled by quoting the
# expansions here instead. Unquoted, one such row word-splits into several
# argv entries and shifts resolution/p1/p2 off the end of the command line,
# which surfaces as run-pwlo.py failing to float() a fragment of the path.
spatial_obj=$(sed -n "${TASK}p" "${spatial_annos}")
out_dir=$(sed -n "${TASK}p" "${outs}")
label=$(sed -n "${TASK}p" "${labels}")

# Echoed with [brackets] before use: a stray CR, a trailing space, or a row
# taken from a misaligned file is invisible in a bare path.
echo "== run-pwlo.s task ${TASK} =="
echo "  spatial_obj : [${spatial_obj}]"
echo "  out_dir     : [${out_dir}]"
echo "  label       : [${label}]"

if [[ ! -f "${spatial_obj}" ]]; then
    echo "run-pwlo.s: spatial object not found: [${spatial_obj}]" >&2
    echo "  Read from row ${TASK} of ${spatial_annos}. Rows must be bare paths," >&2
    echo "  one per line, with no surrounding quotes." >&2
    exit 1
fi

python3 run-pwlo.py "${spatial_obj}" "${out_dir}" "${label}" "${resolution}" "${p1}" "${p2}"
