#!/bin/bash
##SBATCH --mem=50GB
##SBATCH --time=0-2
#SBATCH --error=RunPhenomenalist_%j.err
#SBATCH --out=RunPhenomenalist_%j.out


source phenomenalist-config.txt

module load r/4.1.2

segmentation_file_tmp=$(head -n${SLURM_ARRAY_TASK_ID} ${segmentation_file} | tail -n1)
echo ${segmentation_file_tmp}

failed_markers=$(head -n${SLURM_ARRAY_TASK_ID} ${failed_markers} | tail -n1)
echo ${failed_markers}

nuclear_markers=$(head -n${SLURM_ARRAY_TASK_ID} ${nuclear_markers} | tail -n1)
echo ${nuclear_markers}

out_dir=$(head -n${SLURM_ARRAY_TASK_ID} ${out_dir} | tail -n1)
echo ${out_dir}
mkdir -p ${out_dir}

classifier_labels=$(head -n${SLURM_ARRAY_TASK_ID} ${classifier_labels} | tail -n1)
echo ${classifier_labels}


#names(arguments)=c('segmentation_file','failed.markers','nuclear.markers','HALO','out_dir','clustering_res','classifier_label')

Rscript run-phenomenalist.R ${segmentation_file_tmp} ${failed_markers} ${nuclear_markers} ${HALO} ${out_dir} ${clustering_res} ${classifier_labels} ${max_cells} ${phenotyping_template}
 
