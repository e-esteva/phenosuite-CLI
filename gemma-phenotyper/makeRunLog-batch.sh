#!/bin/bash

file_name=run-log

current_time=$(date "+%Y.%m.%d-%H.%M.%S")

new_fileName=${file_name}.${current_time}.txt

touch ${new_fileName}

source gemma-phenotyper-config.txt

cat gemma-phenotyper-config.txt >> ${new_fileName}
echo 'spe files:' >> ${new_fileName}
cat ${spe_file} >> ${new_fileName}
echo 'out dirs:' >> ${new_fileName}
cat ${out_dir} >> ${new_fileName}
echo 'labels:' >> ${new_fileName}
cat ${labels} >> ${new_fileName}
echo 'cluster cols:' >> ${new_fileName}
cat ${cluster_cols} >> ${new_fileName}
echo 'tissues:' >> ${new_fileName}
cat ${tissues} >> ${new_fileName}

mkdir -p run-logs-batch

mv ${new_fileName} run-logs-batch/
