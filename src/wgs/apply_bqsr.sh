#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BQSR=Null
BAM_FILE=Null
VERSION_LOG=""
OUTPUT=Null
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -bqsr | --bqsr_recal_file | -BQSR )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BQSR="${!OPTARG}"
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
      usage_apply_bqsr
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

[[ ${BQSR} != "Null" ]] || { echo "
ERROR: BQSR (-bqsr <arg>) argument must be provided" && \
 usage_apply_bqsr && exit 1; }
[[ ${BAM_FILE} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_apply_bqsr && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_apply_bqsr && exit 1; }

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/input_data ]] || { MISSING_VOLUMES+=(/data/input_data) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/input_data Read "${BQSR}" || exit 1
python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

    echo "apply_bqsr

Command:
  gatk ApplyBQSR -I /data/bam_files/\"${BAM_FILE}\" \\
    -bqsr /data/input_data/\"${BQSR}\" -O /data/output_data/\"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  gatk:
    version $( get_conda_version gatk )
" > /data/output_data/"${VERSION_LOG}"

fi

gatk ApplyBQSR -I /data/bam_files/"${BAM_FILE}" \
-bqsr /data/input_data/"${BQSR}" -O /data/output_data/"${OUTPUT}"