#! /bin/bash

PILEUP=Null
OUTPUT=Null

while getopts "p:o:h" opt; do
  case ${opt} in
    p )
      PILEUP=${OPTARG}
      ;;
    o )
      OUTPUT=${OPTARG}
      ;;
    h )
      usage_strelka
      exit 0
      ;;
    \? )
      echo "Invalid option: -${OPTARG}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${PILEUP} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_strelka && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_strelka && exit 1; }

if [[ ${REF_GENOME: -3} = ".gz" ]]; then
    INDEX=$(echo ${REF_GENOME} | grep -o '\.' | grep -c '\.')
    NEW_REF="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})".bgz
    gunzip -c /data/ref_index/${REF_GENOME} | bgzip > /data/ref_index/${NEW_REF}
    REF_GENOME=${NEW_REF}
fi

EXIT_CODE=0
NEEDED_FILE=/data/ref_index/${REF_GENOME}.fai
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1

if [[ ! -f ${NEEDED_FILE} ]]; then
    echo "Samtools reference index (${NEEDED_FILE}) is missing. Running samtools faidx"
    samtools faidx /data/ref_index/${REF_GENOME}
fi

varscan somatic ${PILEUP} ${OUTPUT}