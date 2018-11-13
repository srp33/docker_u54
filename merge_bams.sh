#! /bin/bash

BAM_FILES=()
OUTPUT=Null
THREADS=1

while getopts "b:t:o:h" opt; do
  case ${opt} in
    b )
      BAM_FILES+=("/data/bam_files/${OPTARG}")
      ;;
    t )
      THREADS=${OPTARG}
      ;;
    o )
      OUTPUT=${OPTARG}
      ;;
    h )
      usage_merge_bams
      exit 0
      ;;
    \? )
      echo "Invalid option: -${OPTARG}" 1>&2
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
[[ -d /data/results ]] || { MISSING_VOLUMES+=(/data/results) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/results ReadWrite || exit 1


sambamba merge -t ${THREADS} /data/bam_files/${OUTPUT} ${BAM_FILES[@]}