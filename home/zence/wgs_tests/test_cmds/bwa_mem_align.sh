#! /bin/bash

docker run \
  -v /Applications/U54/ref_genome:/data/ref_genome \
  -v /Applications/U54/ref_index:/data/ref_index \
  -v /Applications/U54/SampleData:/data/input_data \
  -v /Applications/U54/output\ data:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  bwa_mem_align \
    -r ucsc.hg19.fasta.gz \
    -s1 101024.1.fastq.gz \
    -s2 101024.2.fastq.gz \
    -o 101024.zence.bam

docker run \
  -v /Applications/U54/ref_genome:/data/ref_genome \
  -v /Applications/U54/ref_index:/data/ref_index \
  -v /Applications/U54/SampleData:/data/input_data \
  -v /Applications/U54/output\ data:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  bwa_mem_align \
    -r ucsc.hg19.fasta.gz \
    -s1 101025.1.fastq.gz \
    -s2 101025.2.fastq.gz \
    -o 101025.zence.bam