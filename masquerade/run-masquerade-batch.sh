#!/bin/bash
#SBATCH --mem=300GB
#SBATCH --partition=cpu_dev
#SBATCH --error=MakeMasks_%j.err
#SBATCH --out=MakeMasks_%j.out

source activate masquerade
module load java/17.0.0

source configFile-batch.txt

image=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${image}")
echo ${image}

spatial_metadata=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${spatial_metadata}")
echo ${spatial_metadata}

outPath=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${outPath}")
echo ${outPath}

relevant_markers=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${relevant_markers}")
echo ${relevant_markers}

mkdir -p run-logs

# Step 1: Generate masks — outputs a .ome.tiff with manual OME-XML
python masquerade_interface.py ${image} ${spatial_metadata} ${outPath} ${relevant_markers} ${adjust_coords} ${compress} ${radius} ${filled} ${num_points} ${preFilter_masks} ${target_size}

# Step 2: Convert large files to pyramidal OME-TIFF via bfconvert if >= 4GB
# Small files (<4GB) are written as .ome.tiff and QuPath reads them directly
out_dir=$(dirname ${outPath})
base=$(basename ${spatial_metadata} .csv)

written_tiff=$(ls -t ${out_dir}/${base}*.ome.tiff 2>/dev/null | head -n1)

if [ -z "${written_tiff}" ]; then
    echo "ERROR: Could not find output .ome.tiff for bfconvert step"
    exit 1
fi

file_size=$(stat -c%s "${written_tiff}")
four_gb=4294967296

if [ "${file_size}" -ge "${four_gb}" ]; then
    pyramidal_out="${written_tiff%.ome.tiff}.pyramidal.ome.tiff"
    echo "File is >= 4GB (${file_size} bytes), converting to pyramidal OME-TIFF: ${pyramidal_out}"

    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

    ${SCRIPT_DIR}/bftools/bfconvert \
        -tilex 512 \
        -tiley 512 \
        -pyramid-resolutions 5 \
        -pyramid-scale 2 \
        -compression LZW \
        "${written_tiff}" \
        "${pyramidal_out}" 2>&1

    if [ $? -eq 0 ]; then
        echo "Conversion successful: ${pyramidal_out}"
        echo "Removing intermediate file: ${written_tiff}"
        rm "${written_tiff}"
        echo "Final output: ${pyramidal_out}"
    else
        echo "ERROR: bfconvert failed — intermediate file retained: ${written_tiff}"
        exit 1
    fi
else
    echo "File is < 4GB (${file_size} bytes), skipping bfconvert"
    echo "Final output: ${written_tiff}"
fi

bash makeRunLog-batch.sh
