#!/bin/bash
#SBATCH --time=0-1
#SBATCH --mem=10GB
#SBATCH --error=RunPhenomenalist-meta_%j.err
#SBATCH --out=RunPhenomenalist-meta_%j.out


source phenomenalist-config.txt 

function phenomenalist {
	local job_id=($(sbatch --export=configfile=$configFile --mem=${module1_mem} --time=${module1_time}  --array=1-${batch_size} --partition=${module1_partition} ${module1_Path}))
	echo ${job_id[3]}
}


echo Initiating phenomenalist at `date`
mod1_job=$(phenomenalist)
echo ${mod1_job}
echo Finished data at `date`
