#! /bin/bash

source usage_functions
source check_functions

set -o errexit

REF_GENOME=Null
KNOWN_SITES=Null
BAM_FILE=Null
VERSION_LOG=""
OUTPUT=Null
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -s | --known_sites )
      check_args "${!OPTARG}" "${!i}" || exit 1
      KNOWN_SITES="${!OPTARG}"
      i=$((i+1))
      ;;
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    -b | --bam_file )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT="${!OPTARG}"
      i=$((i+1))
      ;;
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_base_recalibrator
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

MISSING_VOLUMES=()
EXIT_CODE=0

[[ ${KNOWN_SITES} != "Null" ]] || { echo "
ERROR: KNOWN_SITES (-s <arg>) argument must be provided" && \
 usage_base_recalibrator && exit 1; }
[[ ${REF_GENOME} != "Null" ]] || { echo "
ERROR: REF_GENOME (-r <arg>) argument must be provided" && \
 usage_base_recalibrator && exit 1; }
[[ ${BAM_FILE} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_base_recalibrator && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_base_recalibrator && exit 1; }

[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

#python /check_permissions.py /data/ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1
python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"

INDEX=$(echo "${REF_GENOME}" | grep -o '\.' | grep -c '\.')
if [[ ${REF_GENOME: -${INDEX}} = ".gz" ]]; then
    NEW_REF="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})"
    gunzip -c /data/ref_genome/"${REF_GENOME}" > /tmp/"${NEW_REF}"
    REF_GENOME="${NEW_REF}"
    INDEX=$(echo "${REF_GENOME}" | grep -o '\.' | grep -c '\.')
fi

NEEDED_INDEX=/data/ref_index/"${REF_GENOME}".fai
REF_DICT="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})".dict
NEEDED_DICT=/data/ref_index/"${REF_DICT}"

if [[ ! -f "${NEEDED_INDEX}" ]]; then
    echo "
    Samtools reference index (${NEEDED_INDEX}) is missing. Running samtools faidx
"
    samtools faidx /tmp/"${REF_GENOME}"
    mv /tmp/"${REF_GENOME}.fai" /data/ref_index/"${REF_GENOME}.fai"
fi

ln -s /data/ref_index/"${REF_GENOME}.fai" /tmp/"${REF_GENOME}.fai"

if [[ ! -f ${NEEDED_DICT} ]]; then
    echo "
    Fasta dict file (${NEEDED_DICT}) is missing. Running gatk CreateSequenceDictionary
"
    gatk CreateSequenceDictionary --REF_GENOME /tmp/"${REF_GENOME}"
    mv /tmp/"${REF_DICT}" /data/ref_index/"${REF_DICT}"
fi

ln -s /data/ref_index/"${REF_DICT}" /tmp/"${REF_DICT}"

cd /tmp/
wget -q ${KNOWN_SITES}
cd /data/

INDEX=$(echo ${KNOWN_SITES} | grep -o '/' | grep -c '/')
INDEX=$((INDEX+1))
KNOWN_SITES="$(echo ${KNOWN_SITES} | cut -d '/' -f ${INDEX})"

if [[ ${KNOWN_SITES: -3} = ".gz" ]]; then
    NEW_SITES="$(echo ${KNOWN_SITES} | cut -d '.' -f -3)"
    gunzip -c /tmp/"${KNOWN_SITES}" > /tmp/"${NEW_SITES}"
    KNOWN_SITES="${NEW_SITES}"
    INDEX=$(echo ${KNOWN_SITES} | grep -o '\.' | grep -c '\.')
fi

gatk IndexFeatureFile -F /tmp/${KNOWN_SITES}

if [[ ${VERSION_LOG} != "" ]]; then

    echo "base_recalibrator

Command:
  gatk BaseRecalibrator -R /tmp/\"${REF_GENOME}\" -I /data/bam_files/\"${BAM_FILE}\" \\
    --known-sites /tmp/\"${KNOWN_SITES}\" -O /data/output_data/\"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  gatk:
    version $( get_conda_version gakt )
" > /data/output_data/"${VERSION_LOG}"

fi

gatk BaseRecalibrator -R /tmp/"${REF_GENOME}" -I /data/bam_files/"${BAM_FILE}" \
--known-sites /tmp/"${KNOWN_SITES}" -O /data/output_data/"${OUTPUT}"