#!/bin/bash
#SBATCH --error=pwlo_%j.err
#SBATCH --out=pwlo_%j.out


module load condaenvs/gpu/machinelearning


source config-spatial_dynamics.txt

spatial_obj=$(head -n ${SLURM_ARRAY_TASK_ID} ${spatial_annos} | tail -n1)
echo ${spatial_obj}

out_dir=$(head -n ${SLURM_ARRAY_TASK_ID} ${outs} | tail -n1)
echo ${out_dir}

label=$(head -n ${SLURM_ARRAY_TASK_ID} ${labels} | tail -n1)
echo ${label}

python3 run-pwlo.py ${spatial_obj} ${out_dir} ${label} ${resolution} ${p1} ${p2}
