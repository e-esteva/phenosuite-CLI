#!/bin/bash
#SBATCH --mem=200GB
#SBATCH --time=0-12
#SBATCH --partition=data_mover
#SBATCH --out=copyToSharedDrive_%j.out
#SBATCH --error=copyToSharedDrive_%j.err

mount /mnt/ee699/research/

#destination=/mnt/ee699/research/reizislab/reizislabspace/Craig/CODEX/Processed_QPTIFFs/AMP-output/segmentation/CLAHE-norm/CLAHE/Nuclear+Membrane/out-phenomenalist/AMP_1154_r/

#destination=/mnt/ee699/research/reizislab/reizislabspace/AnnaEichinger/Ergebnisse_Anna/Breast Cancer/Subgated segmentation/

# for path with spaces:
#destination="/mnt/ee699/research/reizislab/reizislabspace/Will/BC-Project/BC_data_for_Eduardo/Gating strategy test/20250212_newstrategy/"
# rsync -avz ${src_} "${destination}"

# -avz to recursively copy contents of folder
rsync -t AMP_1154_r-0.5x.tiff ${destination}

umount /mnt/ee699/research/
