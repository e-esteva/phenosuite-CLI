#!/bin/bash
#SBATCH --mem=10GB
#SBATCH --time=0-1
#SBATCH --error=masquerade_%j.err
#SBATCH --out=masquerade_%j.out


# run batch mask generation:
source configFile-batch.txt

function generate_masks_module {
    local job_id=($(sbatch --export=configfile=$configFile --array=1-${sample_count} --time=0-${run_time} --mem=${memory_} --partition=${masquerade_partition} ${generate_masks_module_Path}))
    echo ${job_id[3]}
}


### MODULE 1
echo Initiating mask generation:
mod1_job=$(generate_masks_module)
echo Finished generating cluster masks at `date`
