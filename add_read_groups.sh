#! /bin/bash

source usage_functions
source check_functions

SAMPLE=Null
GROUP_ID=Null
GROUP_LB=Null
BAM_FILE=Null
VERSION_LOG=""
OUTPUT=Null
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -id | --group_id )
      check_args "${!OPTARG}" "${!i}" || exit 1
      GROUP_ID=${!OPTARG}
      i=$((i+1))
      ;;
    -lb | --library_identifier )
      check_args "${!OPTARG}" "${!i}" || exit 1
      GROUP_LB=${!OPTARG}
      i=$((i+1))
      ;;
    -s | --sample )
      check_args "${!OPTARG}" "${!i}" || exit 1
      SAMPLE="${!OPTARG}"
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
      usage_add_read_groups
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

[[ ${GROUP_ID} != "Null" ]] || { echo "
ERROR: ID (-id <arg>) argument must be provided" && \
 usage_add_read_groups && exit 1; }
[[ ${GROUP_LB} != "Null" ]] || { echo "
ERROR: LB (-lb <arg>) argument must be provided" && \
 usage_add_read_groups && exit 1; }
[[ ${SAMPLE} != "Null" ]] || { echo "
ERROR: SAMPLE (-s <arg>) argument must be provided" && \
 usage_add_read_groups && exit 1; }
[[ ${BAM_FILE} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_add_read_groups && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_add_read_groups && exit 1; }

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

    echo "add_read_groups

Command:
  samtools addreplacerg -r ID:\"${GROUP_ID}\" -r LB:\"${GROUP_LB}\" -r SM:\"${SAMPLE}\" \\
    -o /data/output_data/\"${OUTPUT}\" /data/bam_files/\"${BAM_FILE}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  samtools:
    version $( get_conda_version samtools )
" > /data/output_data/"${VERSION_LOG}"

fi

samtools addreplacerg -r ID:"${GROUP_ID}" -r LB:"${GROUP_LB}" -r SM:"${SAMPLE}" \
-o /data/output_data/"${OUTPUT}" /data/bam_files/"${BAM_FILE}"