#!/bin/bash
#SBATCH --mem=10GB
#SBATCH --time=0-1
#SBATCH --error=spatial_dynamics_meta_%j.err
#SBATCH --out=spatial_dynamics_meta_%j.out


# run batch PWLO:
source config-spatial_dynamics.txt

function run_pwlo_module {
    local job_id=($(sbatch --export=configfile=$configFile --array=1-${sample_count} --time=0-${run_time} --mem=${memory_} --partition=${spatial_dynamics_meta_partition} ${pwlo_module_Path}))
    echo ${job_id[3]}
}

function run_spatial_circuit_module {
    local job_id=($(sbatch --export=configfile=$configFile --array=1-${sample_count} --time=0-${run_time} --mem=${memory_} --partition=${spatial_dynamics_meta_partition} ${spatial_circuit_module_Path}))
    echo ${job_id[3]}
}

if [[ ${module} -eq 0 ]]
then
	### MODULE 1
	echo Initiating pairwise log-odds calculation:
	mod1_job=$(run_pwlo_module)
	echo Finished at `date`
else
	### MODULE 2
	echo Initiating multicellular spatial circuit enrichment calculation:
	mod1_job=$(run_spatial_circuit_module)
	echo Finished at `date`

fi
