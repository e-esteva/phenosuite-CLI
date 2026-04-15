#!/bin/bash

file_name=run-log

current_time=$(date "+%Y.%m.%d-%H.%M.%S")

new_fileName=${file_name}.${current_time}.txt

touch ${new_fileName}

source configFile-batch.txt

cat configFile-batch.txt >> ${new_fileName}

cat ${spatial_metadata} >> ${new_fileName}
cat ${relevant_markers} >> ${new_fileName}
cat ${outDir} >> ${new_fileName}
cat ${outPath} >> ${new_fileName}
cat ${image} >> ${new_fileName}

mkdir -p run-logs-batch

mv ${new_fileName} run-logs-batch/

