#! /bin/bash

source usage_functions
source check_functions

set -o errexit

MAX_DISTANCE=1000
MIN_AGREEMENT=2
ACCOUNT_TYPE=1
ACCOUNT_STRANDS=1
BASE_ON_SIZE=0
MIN_SIZE=30
VCF_ARG=Null
OUTPUT=Null
SAMPLE=""
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -v | --input_vcfs )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VCF_ARG="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT="${!OPTARG}"
      i=$((i+1))
      ;;
    -max | --max_distance )
      check_args "${!OPTARG}" "${!i}" || exit 1
      MAX_DISTANCE="${!OPTARG}"
      i=$((i+1))
      ;;
    -s | --sample )
      check_args "${!OPTARG}" "${!i}" || exit 1
      SAMPLE="${!OPTARG}"
      i=$((i+1))
      ;;
    --min_agreement )
      check_args "${!OPTARG}" "${!i}" || exit 1
      MIN_AGREEMENT="${!OPTARG}"
      i=$((i+1))
      ;;
    --min_size )
      check_args "${!OPTARG}" "${!i}" || exit 1
      MIN_SIZE="${!OPTARG}"
      i=$((i+1))
      ;;
    --ignore_type )
      ACCOUNT_TYPE=0
      ;;
    --ignore_strand )
      ACCOUNT_STRANDS=0
      ;;
    --estimate_by_size )
      BASE_ON_SIZE=1
      ;;
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_survivor
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${VCF_ARG} != "Null" ]] || { echo "
ERROR: VCF FILE (-i <arg>) argument must be provided" && \
 usage_survivor && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_survivor && exit 1; }

if [[ -d /data/input_data/"${VCF_ARG}" ]]; then
  for vcf_file in /data/input_data/"${VCF_ARG}"/"${SAMPLE}"*.vcf; do
    ln -s /data/input_data/"${VCF_ARG}"/"${vcf_file##*/}" ./"${vcf_file##*/}"
  done
  ls ${SAMPLE}*.vcf > sample_files
  #chmod 755 sample_files
  VCF_FILE=sample_files
elif [[ -f /data/input/data/"${VCF_ARG}" ]]; then
  echo "ERROR: File with VCF names not currently supported. Please give directory containing \
VCF files as input. If necessary, give sample (or common beginning) of desired VCF files."
  exit 1
else
  echo "ERROR: ${VCF_ARG} not found!!"
  exit 1
fi


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

    echo "run_survivor

Command:
  SURVIVOR merge \
    \"${VCF_FILE}\" \
    ${MAX_DISTANCE} \
    ${MIN_AGREEMENT} \
    ${ACCOUNT_TYPE} \
    ${ACCOUNT_STRANDS} \
    ${BASE_ON_SIZE} \
    ${MIN_SIZE} \
    /data/output_data/\"${OUTPUT}\"

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

source activate py2.7

survivor merge \
  "${VCF_FILE}" \
  ${MAX_DISTANCE} \
  ${MIN_AGREEMENT} \
  ${ACCOUNT_TYPE} \
  ${ACCOUNT_STRANDS} \
  ${BASE_ON_SIZE} \
  ${MIN_SIZE} \
  /data/output_data/"${OUTPUT}"