#! /bin/bash

source usage_functions
source check_functions

set -o errexit

TUMOR=Null
NORMAL=Null
REF_GENOME=Null
OVERWRITE=0
VERSION_LOG=""
INDEL=""
CALL_REGIONS=""
RUN_DIR="/data/output_data/StrelkaSomaticWorkflow"
RUN_DIR_ARG="--runDir=${RUN_DIR}"
ARGNUM=$#
EXOME=""

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -t | --tumorBam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      TUMOR="${!OPTARG}"
      i=$((i+1))
      ;;
    -n | --normalBam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      NORMAL="${!OPTARG}"
      i=$((i+1))
      ;;
    -r | --referenceFasta )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    -i | --indelCandidates )
      check_args "${!OPTARG}" "${!i}" || exit 1
      INDEL="--indelCandidates /data/input_data/${!OPTARG}"
      i=$((i+1))
      ;;
    -c | --callRegions )
      check_args "${!OPTARG}" "${!i}" || exit 1
      CALL_REGIONS="--callRegions /data/input_data/${!OPTARG}"
      i=$((i+1))
      ;;
    -d | --runDir )
      check_args "${!OPTARG}" "${!i}" || exit 1
      RUN_DIR=/data/output_data/"${!OPTARG}"
      RUN_DIR_ARG="--runDir=${RUN_DIR}"
      i=$((i+1))
      ;;
    -e | --exome )\
      EXOME="--exome"
      ;;
    -O | --overwrite_runDir )
      OVERWRITE=1
      ;;
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
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

[[ "${TUMOR}" != "Null" ]] || { echo "
ERROR: TUMOR BAM FILE (-t <arg>) argument must be provided" && \
 usage_strelka && exit 1; }
[[ "${NORMAL}" != "Null" ]] || { echo "
ERROR: NORMAL BAM FILE (-n <arg>) argument must be provided" && \
 usage_strelka && exit 1; }
[[ "${REF_GENOME}" != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_strelka && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

[[ -f /data/bam_files/"${TUMOR}.bai" ]] || { echo "
ERROR: /data/bam_files/\"${TUMOR}.bai\" does not exist. Please run index_bam on \
/data/bam_files/\"${TUMOR}\" in order to use this command.
" && usage_strelka && exit 1; }
[[ -f /data/bam_files/"${NORMAL}.bai" ]] || { echo "
ERROR: /data/bam_files/\"${NORMAL}.bai\" does not exist. Please run index_bam on \
/data/bam_files/\"${NORMAL}\" in order to use this command.
" && usage_strelka && exit 1; }

#mkdir /temp
ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"

INDEX=$(echo "${REF_GENOME}" | grep -o '\.' | grep -c '\.')
if [[ "${REF_GENOME: -${INDEX}}" = ".gz" ]]; then
    NEW_REF="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})"
    gunzip -c /data/ref_genome/"${REF_GENOME}" > /tmp/"${NEW_REF}"
    REF_GENOME="${NEW_REF}"
fi

NEEDED_FILE=/data/ref_index/"${REF_GENOME}".fai

if [[ ! -f "${NEEDED_FILE}" ]]; then
    echo "
    Samtools reference index (${NEEDED_FILE}) is missing. Running samtools faidx
"
    samtools faidx /tmp/"${REF_GENOME}"
    mv /tmp/"${REF_GENOME}.fai" /data/ref_index/"${REF_GENOME}.fai"
fi

ln -s /data/ref_index/"${REF_GENOME}.fai" /tmp/"${REF_GENOME}.fai"

if [[ ${VERSION_LOG} != "" ]]; then

    echo "call_somatic_variants_strelka

Commands:
  python2.7 /miniconda/share/strelka-2.9.10-0/system_commands/configureStrelkaSomaticWorkflow.py \\
    --referenceFasta=/tmp/\"${REF_GENOME}\" --tumorBam=\"/data/bam_files/${TUMOR}\" \\
    --normalBam=\"/data/bam_files/${NORMAL}\" ${INDEL} ${CALL_REGIONS} \"${RUN_DIR_ARG}\" \\
    \"${EXOME}\"

  python2.7 \"${RUN_DIR}\"/runWorkflow.py --mode=local

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python2_version )

  samtools:
    version $( get_conda_version samtools )

  strelka:
    version 2.9.10-0
" > /data/output_data/"${VERSION_LOG}"

fi

[[ ${OVERWRITE} = 0 ]] || rm -f "${RUN_DIR_ARG##*=}/runWorkflow.py"

#source activate py2.7

python2.7 /miniconda/envs/py2.7/share/strelka-2.9.10-0/bin/configureStrelkaSomaticWorkflow.py \
    --referenceFasta=/tmp/"${REF_GENOME}" --tumorBam="/data/bam_files/${TUMOR}" \
    --normalBam="/data/bam_files/${NORMAL}" ${INDEL} ${CALL_REGIONS} "${RUN_DIR_ARG}" \
    ${EXOME}

python2.7 "${RUN_DIR}"/runWorkflow.py --mode=local
