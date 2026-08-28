#!/bin/bash
#SBATCH --mem=200GB
#SBATCH --time=0-2
#SBATCH --error=spatial_circuit_%j.err
#SBATCH --out=spatial_circuit_%j.out

# Activate the shipped conda env. Provision once with:
#   conda env create -f environment.yml
source activate spatial-dynamics

source config-spatial_dynamics.txt

TASK="${SLURM_ARRAY_TASK_ID}"

# Same contract as run-pwlo.s: rows are bare paths, expansions are quoted
# here so paths containing spaces survive as a single argument.
spatial_obj=$(sed -n "${TASK}p" "${spatial_annos}")
out_dir=$(sed -n "${TASK}p" "${outs}")
label=$(sed -n "${TASK}p" "${labels}")

echo "== run-spatial_circuit-enrichment.s task ${TASK} =="
echo "  spatial_obj : [${spatial_obj}]"
echo "  out_dir     : [${out_dir}]"
echo "  label       : [${label}]"
echo "  circuit     : [${circuit}]"

if [[ ! -f "${spatial_obj}" ]]; then
    echo "run-spatial_circuit-enrichment.s: spatial object not found: [${spatial_obj}]" >&2
    echo "  Read from row ${TASK} of ${spatial_annos}. Rows must be bare paths," >&2
    echo "  one per line, with no surrounding quotes." >&2
    exit 1
fi

python3 run-spatial_circuit-enrichment.py "${spatial_obj}" "${out_dir}" "${label}" "${circuit}"
