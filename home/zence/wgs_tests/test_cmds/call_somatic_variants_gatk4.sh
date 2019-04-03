#! /bin/bash

docker run \
  -v /Applications/U54/bam_files:/data/bam_files \
  -v /Applications/U54/ref_genome:/data/ref_genome \
  -v /Applications/U54/ref_index:/data/ref_index \
  -v /Applications/U54/output\ data:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  call_somatic_variants_gatk4 \
    -t 101025.zence.bam \
    -ts 101025 \
    -n 101024.zence.bam \
    -ns 101024 \
    -r ucsc.hg19.fasta.gz \
    -o 101024.101025.zence.bam