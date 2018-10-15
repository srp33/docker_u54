#! /bin/bash

REF_GENOME=$1
BWA_THREADS=$2
SAMPLES=()
FILE_EXTENSIONS=()
NEEDED_FILES=(amb ann bwt pac sa)
REF_INDEXED=0

for filename in /data/ref_index/*; do
	echo "${filename##*.}"
    FILE_EXTENSIONS+=($(echo "${filename##*.}"))
    echo ${FILE_EXTENSIONS[@]}
done

for NEEDED_EXT in ${NEEDED_FILES[@]}; do
    [[ " ${FILE_EXTENSIONS[@]} " =~ " ${NEEDED_EXT} " ]] || REF_INDEXED=1
done

echo ${REF_INDEXED}

if [[ REF_INDEXED == 1 ]]; then

    read -p "The reference does not contain the proper index files. Run bwa index? [y/n]" -n 1 -s to_index

    if [ "$to_index" -eq "y" ]; then
        bwa index /data/ref_index/${REF_GENOME}
    else
        echo "Index files for reference genome are required to run bwa mem"
        exit 1
    fi
fi

for i in /data/sample_data/*; do
    SAMPLES+=($(echo ${i} | cut -d '/' -f 4 |  cut -d '.' -f 1))
done

SAMPLES=($(echo ${SAMPLES[@]} | tr ' ' '\n' | sort -u | tr '\n' ' '))

for sample in ${SAMPLES[@]}; do
    bwa mem -t ${BWA_THREADS} /data/ref_index/${REF_GENOME} \
        /data/sample_data/${sample}.1.fastq.gz /data/sample_data/${sample}.2.fastq.gz | \
        samtools view -S -b > /data/results/${sample}.bam
done
