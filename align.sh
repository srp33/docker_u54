#! /bin/bash

REF_GENOME=Null
THREADS=1
SAMPLE=Null

while getopts "t:r:s:" opt; do
  case ${opt} in
    t )
      THREADS=${OPTARG}
      ;;
    r )
      REF_GENOME=${OPTARG}
      ;;
    s )
      SAMPLE=${OPTARG}
      ;;
    \? )
      echo "Invalid option: -${OPTARG}" 1>&2
      exit 1
      ;;
  esac
done


REF_INDEX_FILES=()
NEEDED_FILES=(${REF_GENOME}.amb ${REF_GENOME}.ann ${REF_GENOME}.bwt ${REF_GENOME}.pac ${REF_GENOME}.sa)
REF_INDEXED=0
MISSING_VOLUMES=()
EXIT_CODE=0

# Since the default value for REF_GENOME and SAMPLE is "Null", we can check if those have been changed

[[ ${REF_GENOME} != "Null" ]] || { echo "ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && exit 1; }
[[ ${SAMPLE} != "Null" ]] || { echo "ERROR: SAMPLE (-s <arg>) argument must be provided" && exit 1; }

# Checks for the necessary directories which are only created by volumes

[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/sample_data ]] || { MISSING_VOLUMES+=(/data/sample_data) && EXIT_CODE=1; }
[[ -d /data/results ]] || { MISSING_VOLUMES+=(/data/results) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

# Check permissions of each directory

python /check_permissions.py /data/ref_index Read ${REF_GENOME} || exit 1
python /check_permissions.py /data/sample_data Read ${SAMPLE}.1.fastq.gz || exit 1
python /check_permissions.py /data/results ReadWrite || exit 1

# Check for necessary index files in ref_index directory
#   If one of the files is missing, bwa index will be run

for filename in /data/ref_index/*; do
    REF_INDEX_FILES+=($(echo "${filename##*/}"))
done

for NEEDED_FILE in ${NEEDED_FILES[@]}; do
    [[ " ${REF_INDEX_FILES[@]} " =~ " ${NEEDED_FILE} " ]] || REF_INDEXED=1
done


if [[ ${REF_INDEXED} == 1 ]]; then

    echo "The reference does not contain the proper index files. Running bwa index"

    python check_permissions.py /data/ref_index ReadWrite || \
       { echo "Please ensure you are passing in directory and not just a file volume" && exit 1; }

    bwa index -t ${THREADS} /data/ref_index/${REF_GENOME}
fi

bwa mem -t ${THREADS} /data/ref_index/${REF_GENOME} \
    /data/sample_data/${SAMPLE}.1.fastq.gz /data/sample_data/${SAMPLE}.2.fastq.gz | \
    samtools view -@ ${THREADS} -S -b > /data/results/${SAMPLE}.bam

