#! /bin/bash

REF_GENOME=$1
BWA_THREADS=$2
SAMPLES=()

pwd

for i in /data/SampleData/*; do
        echo ${i}
        SAMPLES+=($(echo ${i} | cut -d '/' -f 4 |  cut -d '.' -f 1))
done
echo ${SAMPLES[@]}
SAMPLES=($(echo ${SAMPLES[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '))

for sample in ${SAMPLES[@]}; do
        echo ${sample}
        bwa mem -t ${BWA_THREADS} /data/ref_index/${REF_GENOME} \
        /data/SampleData/${sample}.1.fastq.gz /data/SampleData/${sample}.2.fastq.gz | \
            samtools view -S -b > /data/results/${sample}.bam
done
