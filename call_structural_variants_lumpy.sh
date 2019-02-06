#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BAM_FILE=Null
DISC_FILE=Null
SPLIT_FILE=Null
OUTPUT=Null
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -b | --bam_file )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -s | --split_reads )
      check_args "${!OPTARG}" "${!i}" || exit 1
      SPLIT_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -d | --discordant_reads )
      check_args "${!OPTARG}" "${!i}" || exit 1
      DISC_FILE="${!OPTARG}"
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
      usage_manta
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ "${BAM_FILE}" != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_manta && exit 1; }
[[ "${DISC_FILE}" != "Null" ]] || { echo "
ERROR: DISCORDANT READS SAM FILE (-d <arg>) argument must be provided" && \
 usage_manta && exit 1; }
[[ "${SPLIT_FILE}" != "Null" ]] || { echo "
ERROR: SPLIT READS SAM FILE (-s <arg>) argument must be provided" && \
 usage_manta && exit 1; }
[[ "${OUTPUT}" != "Null" ]] || { echo "
ERROR: OUTPUT (-s <arg>) argument must be provided" && \
 usage_manta && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/input_data ]] || { MISSING_VOLUMES+=(/data/input_data) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/input_data Read "${SPLIT_FILE}" || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

    echo "call_somatic_variants_strelka

Commands:
  lumpyexpress -B /data/bam_files/\"${BAM_FILE}\" -S /data/input_data/\"${SPLIT_FILE}\" \\
    -D /data/input_data/\"${DISC_FILE}\" -o /data/output_data/\"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python2_version )

  lumpy:
    version $( get_conda_version lumpy )

  strelka:
    version 2.9.10-0
" > /data/output_data/"${VERSION_LOG}"

fi

source activate py2.7

lumpyexpress -B /data/bam_files/"${BAM_FILE}" -S /data/input_data/"${SPLIT_FILE}" \
    -D /data/input_data/"${DISC_FILE}" -o /data/output_data/"${OUTPUT}"
