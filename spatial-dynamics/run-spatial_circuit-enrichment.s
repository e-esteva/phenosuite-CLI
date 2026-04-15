#!/bin/bash
#SBATCH --mem=200GB
#SBATCH --time=0-2
#SBATCH --error=spatial_circuit_%j.err
#SBATCH --out=spatial_circuit_%j.out

module load condaenvs/gpu/machinelearning


source config-spatial_dynamics.txt

python3 run-spatial_circuit-enrichment.py ${spatial_obj} ${out_dir} ${label} ${circuit}
