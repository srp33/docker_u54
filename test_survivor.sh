#! /bin/bash

source activate py2.7

echo "Now running in $CONDA_DEFAULT_ENV"

which SURVIVOR

cd /data/input_data/piped_vcfs/101024

ls /data/input_data/piped_vcfs/101024*.vcf > sample_files

SURVIVOR merge sample_files 1000 2 1 1 0 30 testing.surv.vcf