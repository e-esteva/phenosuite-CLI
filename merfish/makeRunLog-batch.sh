#!/bin/bash
#
# makeRunLog-batch.sh — snapshot the config + batch-input lists into a timestamped
# run log under run-logs-batch/, so each submission is reproducible after the fact.

file_name=run-log
current_time=$(date "+%Y.%m.%d-%H.%M.%S")
new_fileName=${file_name}.${current_time}.txt

touch "${new_fileName}"

source merfish-config.txt

cat merfish-config.txt >> "${new_fileName}"
echo ''                >> "${new_fileName}"
echo 'expression files:' >> "${new_fileName}"; cat "${expression_file}" >> "${new_fileName}"
echo 'metadata files:'   >> "${new_fileName}"; cat "${metadata_file}"   >> "${new_fileName}"
echo 'out dirs:'         >> "${new_fileName}"; cat "${out_dir}"         >> "${new_fileName}"
if [[ -n "${sample_id:-}" && -f "${sample_id}" ]]; then
    echo 'sample ids:'   >> "${new_fileName}"; cat "${sample_id}"       >> "${new_fileName}"
fi

mkdir -p run-logs-batch
mv "${new_fileName}" run-logs-batch/
echo "Wrote run-logs-batch/${new_fileName}"
