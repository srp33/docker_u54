#! /bin/bash

docker run \
  -v /Application/U54/bam_files:/data/bam_files \
  -v /Application/U54/output\ data/:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  add_read_groups \
    -b 101024.zence.bam \
    -id 101024 \
    -lb zence \
    -s 101024 \
    -o 101024.zence.bam
