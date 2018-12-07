#! /bin/bash

source usage_functions
source check_functions

BAM_FILE=Null
OUTPUT=Null
VERSION_LOG=""
THREADS=1
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -b | --bam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE=${!OPTARG}
      i=$((i+1))
      ;;
    -t | --nthreads )
      check_args "${!OPTARG}" "${!i}" || exit 1
      THREADS=${!OPTARG}
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT=${!OPTARG}
      i=$((i+1))
      ;;
    --version )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_sort_bam
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
 usage_sort_bam && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_sort_bam && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

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
  sambamba sort -t ${THREADS} -o /data/output_data/\"${OUTPUT}\" /data/bam_files/\"${BAM_FILE}\"

Date run: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  samtools:
    version $( get_conda_version samtools )
" > /data/output_data/"${VERSION_LOG}"

fi

sambamba sort -t ${THREADS} -o /data/output_data/"${OUTPUT}" /data/bam_files/"${BAM_FILE}"