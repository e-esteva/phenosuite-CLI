#!/bin/bash
# makeRunLog-batch.sh
# Creates a timestamped run log capturing the current config + batch inputs.

file_name=run-log
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
new_fileName=${file_name}.${current_time}.txt

touch "${new_fileName}"

source pcf-config.txt

cat pcf-config.txt        >> "${new_fileName}"
echo 'vectra_files:'      >> "${new_fileName}"
cat "${vectra_files}"     >> "${new_fileName}"
echo 'out_dirs:'          >> "${new_fileName}"
cat "${out_dirs}"         >> "${new_fileName}"
echo 'labels:'            >> "${new_fileName}"
cat "${labels}"           >> "${new_fileName}"
echo 'celltypes:'         >> "${new_fileName}"
cat "${celltypes}"        >> "${new_fileName}"
echo 'ref_celltypes:'     >> "${new_fileName}"
cat "${ref_celltypes}"    >> "${new_fileName}"

mkdir -p run-logs-batch
mv "${new_fileName}" run-logs-batch/
