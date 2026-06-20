#!/bin/bash
# makeRunLog-batch.sh
# Creates a timestamped run log capturing the current config + batch inputs.

file_name=run-log
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
new_fileName=${file_name}.${current_time}.txt

touch "${new_fileName}"

source neighborhood_analysis-config.txt

cat neighborhood_analysis-config.txt  >> "${new_fileName}"
echo 'rds_files:'                     >> "${new_fileName}"
cat "${rds_files}"                    >> "${new_fileName}"
echo 'celltype_cols:'                 >> "${new_fileName}"
cat "${celltype_cols}"                >> "${new_fileName}"
echo 'out_dirs:'                      >> "${new_fileName}"
cat "${out_dirs}"                     >> "${new_fileName}"
echo 'labels:'                        >> "${new_fileName}"
cat "${labels}"                       >> "${new_fileName}"
echo 'condition_maps:'                >> "${new_fileName}"
cat "${condition_maps}"               >> "${new_fileName}"

mkdir -p run-logs-batch
mv "${new_fileName}" run-logs-batch/
