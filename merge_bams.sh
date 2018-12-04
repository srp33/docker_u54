#! /bin/bash

source usage_functions
source check_for_args

BAM_FILES=()
OUTPUT=Null
THREADS=1
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -b | --bam )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILES+=("/data/bam_files/${!OPTARG}")
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
    -h | --help )
      usage_merge_bams
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${#BAM_FILES[@]} -gt 1 ]] || { echo "
ERROR: TWO OR MORE BAM FILE (-b <arg>) arguments must be provided" && \
 usage_merge_bams && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) arguments must be provided" && \
 usage_merge_bams && exit 1; }

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


sambamba merge -t ${THREADS} /data/output_data/"${OUTPUT}" ${BAM_FILES[@]}