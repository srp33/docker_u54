#! /bin/bash

docker run \
  -v /Applications/U54/bam_files:/data/bam_files \
  -v /Applications/U54/output\ data:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  call_structural_variants_lumpy \
    -b 101024.zence.bam \
    -s 101024.zence.split.bam \
    -d 101024.zence.disc.bam \
    -o 101024.zence.lumpy.vcf
