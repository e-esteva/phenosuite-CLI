#!/bin/bash
#SBATCH --mem=200GB
#SBATCH --time=0-2
#SBATCH --error=spatial_circuit_%j.err
#SBATCH --out=spatial_circuit_%j.out

# Activate the shipped conda env. Provision once with:
#   conda env create -f environment.yml
source activate spatial-dynamics

source config-spatial_dynamics.txt

python3 run-spatial_circuit-enrichment.py ${spatial_obj} ${out_dir} ${label} ${circuit}
