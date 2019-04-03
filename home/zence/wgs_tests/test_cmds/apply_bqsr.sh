#! /bin/bash

docker run \
  -v /Applications/U54/bam_files:/data/bam_files \
  -v /Applications/U54/in_use:/data/input_data \
  -v /Applications/U54/output\ data:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  apply_bqsr \
    -b 101024.zence.bam \
    -bqsr recal_data.table \
    -o 101024.zence.bam
