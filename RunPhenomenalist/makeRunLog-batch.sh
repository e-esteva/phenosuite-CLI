#!/bin/bash

file_name=run-log

current_time=$(date "+%Y.%m.%d-%H.%M.%S")

new_fileName=${file_name}.${current_time}.txt

touch ${new_fileName}

source phenomenalist-config.txt

cat phenomenalist-config.txt >> ${new_fileName}
echo 'segmentation files:' >> ${new_fileName}
cat ${segmentation_file} >> ${new_fileName}
echo 'failed markers:' >> ${new_fileName}
cat ${failed_markers} >> ${new_fileName}
echo 'nuclear markers:' >> ${new_fileName}
cat ${nuclear_markers} >> ${new_fileName}
echo 'out dirs:' >> ${new_fileName}
cat ${out_dir} >> ${new_fileName}
echo 'classifier labels:' >> ${new_fileName}
cat ${classifier_labels} >> ${new_fileName}

mkdir -p run-logs-batch

mv ${new_fileName} run-logs-batch/

