#!/bin/bash
# makeRunLog-batch.sh
# Creates a timestamped run log capturing the current config + batch inputs.

file_name=run-log
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
new_fileName=${file_name}.${current_time}.txt

touch "${new_fileName}"

source segmentation-config.txt

cat segmentation-config.txt         >> "${new_fileName}"
echo 'image:'                       >> "${new_fileName}"
cat "${image}"                      >> "${new_fileName}"
echo 'out_dirs:'                    >> "${new_fileName}"
cat "${out_dirs}"                   >> "${new_fileName}"
echo 'labels:'                      >> "${new_fileName}"
cat "${labels}"                     >> "${new_fileName}"
echo 'methods:'                     >> "${new_fileName}"
cat "${methods}"                    >> "${new_fileName}"
echo 'nuclear_channels:'            >> "${new_fileName}"
cat "${nuclear_channels}"           >> "${new_fileName}"
echo 'membrane_channels:'           >> "${new_fileName}"
cat "${membrane_channels}"          >> "${new_fileName}"

mkdir -p run-logs-batch
mv "${new_fileName}" run-logs-batch/
