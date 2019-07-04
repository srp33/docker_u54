#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BAM_FILE=Null
VCF_FILE=Null
JSON_FILE=""
OUTPUT=Null
VERSION_LOG=""
THREADS=1
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -b | --bam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -i | --input_vcf )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VCF_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -l | --lib_info )
      check_args "${!OPTARG}" "${!i}" || exit 1
      JSON_FILE="-l /data/input_data/${!OPTARG}"
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
      usage_svtype
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${BAM_FILE} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_svtype && exit 1; }
[[ ${VCF_FILE} != "Null" ]] || { echo "
ERROR: VCF FILE (-i <arg>) argument must be provided" && \
 usage_svtype && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_svtype && exit 1; }

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
python /check_permissions.py /data/input_data Read "${VCF_FILE}" || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

    echo "svtype_vcfs

Command:
  svtyper \
    -i /data/input_data/\"${VCF_FILE}\" \
    -b /data/bam_files/\"${BAM_FILE}\" \
    ${JSON_FILE} \
    > /data/output_data/\"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  svtyper:
    version $( get_conda_version svtyper )
" > /data/output_data/"${VERSION_LOG}"

fi

echo "Running svtyper on ${VCF_FILE}..."

source activate py2.7

svtyper \
    -i /data/input_data/"${VCF_FILE}" \
    -B /data/bam_files/"${BAM_FILE}" \
    ${JSON_FILE} \
    > /data/output_data/"${OUTPUT}"

echo "Done!"
