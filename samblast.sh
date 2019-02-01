#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BAM_FILE=Null
OUTPUT="/dev/null"
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -b | --bam_file )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT="/data/output_data/${!OPTARG}"
      i=$((i+1))
      ;;
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_align
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

[[ ${BAM_FILE} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_samblast && exit 1; }

SAMPLE="${BAM_FILE%%.*}"

# Checks for the necessary directories which are only created by volumes

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

# Check permissions of each directory

python /check_permissions.py /data/bam_files Read "${BAM_FILE}" || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1


if [[ ${VERSION_LOG} != "" ]]; then

    echo "samblast

Command:
  samtools view -h /data/bam_files/\"${BAM_FILE}\" | samblaster -a -e \\
-d /data/output_data/\"${SAMPLE}.disc.sam\" -s /data/output_data/\"${SAMPLE}.split.sam\" \\
-o \"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  samblaster:
    version $( get_conda_version samblaster )

  samtools:
    version $( get_conda_version samtools )
" > /data/output_data/"${VERSION_LOG}"

fi

samtools view -h /data/bam_files/"${BAM_FILE}" | samblaster --addMateTags --ignoreUnmated \
| samblaster -a -e --ignoreUnmated \
-d /data/output_data/"${SAMPLE}.disc.sam" -s /data/output_data/"${SAMPLE}.split.sam" \
-o "${OUTPUT}"