#! /bin/bash

source usage_functions
source check_functions

TUMOR=Null
NORMAL=Null
REF_GENOME=Null
INDEL=""
CALL_REGIONS=""
RUN_DIR="/data/output_data/StrelkaSomaticWorkflow"
RUN_DIR_ARG="--runDir ${RUN_DIR}"
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -t | --tumorBam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      TUMOR=${!OPTARG}
      i=$((i+1))
      ;;
    -n | --normalBam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      NORMAL=${!OPTARG}
      i=$((i+1))
      ;;
    -r | --referenceFasta )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME=${!OPTARG}
      i=$((i+1))
      ;;
    -i | --indelCandidates )
      check_args "${!OPTARG}" "${!i}" || exit 1
      INDEL="--indelCandidates ${!OPTARG}"
      i=$((i+1))
      ;;
    -c | --callRegions )
      check_args "${!OPTARG}" "${!i}" || exit 1
      CALL_REGIONS="--callRegions ${!OPTARG}"
      i=$((i+1))
      ;;
    -d | --runDir )
      check_args "${!OPTARG}" "${!i}" || exit 1
      RUN_DIR=/data/output_data/${!OPTARG}
      RUN_DIR_ARG="--runDir ${RUN_DIR}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_strelka
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${TUMOR} != "Null" ]] || { echo "
ERROR: TUMOR BAM FILE (-t <arg>) argument must be provided" && \
 usage_strelka && exit 1; }
[[ ${NORMAL} != "Null" ]] || { echo "
ERROR: NORMAL BAM FILE (-n <arg>) argument must be provided" && \
 usage_strelka && exit 1; }
[[ ${REF_GENOME} != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_strelka && exit 1; }


EXIT_CODE=0
NEEDED_FILE=/data/ref_genome/${REF_GENOME}.fai
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/ref_genome ReadWrite || exit 1

if [[ ! -f ${NEEDED_FILE} ]]; then
    echo "
    Samtools reference index (${NEEDED_FILE}) is missing. Running samtools faidx
"
    INDEX=$(echo ${REF_GENOME} | grep -o '\.' | grep -c '\.')
    if [[ ${REF_GENOME: -${INDEX}} = ".gz" ]]; then
        mkdir /data/temp_ref
        NEW_REF="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})"
        gunzip -c /data/ref_genome/${REF_GENOME} > /data/temp_ref/${NEW_REF}
        REF_LOCATION=/data/temp_ref/${NEW_REF}
    fi
    samtools faidx ${REF_LOCATION}
fi


python2.7 /opt/miniconda/share/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py \
--referenceFasta="${REF_LOCATION}" --tumorBam="/data/bam_files/${TUMOR}" \
--normalBam="/data/bam_files/${NORMAL}" "${INDEL}" "${CALL_REGIONS}" "${RUN_DIR_ARG}"

python2.7 "${RUN_DIR}"/runWorkflow.py --mode=local
