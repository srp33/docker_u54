#! /bin/bash

REF_GENOME=$1
BWA_THREADS=$2
SAMPLES=()

#read -p "The reference does not contain the proper index files. Run bwa index? [y/n]" -n 1 -s to_index

#if [ "$to_index" -eq "y" ]; then
#    bwa index /data/ref_index/${REF_GENOME}
#else
#    echo "Index files for reference genome are required to run bwa"
#    exit 1
#fi

for i in /data/SampleData/*; do
        SAMPLES+=($(echo ${i} | cut -d '/' -f 4 |  cut -d '.' -f 1))
done

SAMPLES=($(echo ${SAMPLES[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '))

for sample in ${SAMPLES[@]}; do
        bwa mem -t ${BWA_THREADS} /data/ref_index/${REF_GENOME} \
        /data/SampleData/${sample}.1.fastq.gz /data/SampleData/${sample}.2.fastq.gz | \
            samtools view -S -b > /data/results/${sample}.bam
done
