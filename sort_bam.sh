#! /bin/bash

BAM_FILE=Null
THREADS=1

while getopts "b:t:h" opt; do
  case ${opt} in
    b )
      BAM_FILE=${OPTARG}
      ;;
    t )
      THREADS=${OPTARG}
      ;;
    h )
      usage_sort_bam
      exit 0
      ;;
    \? )
      echo "Invalid option: -${OPTARG}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${BAM_FILE} != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_sort_bam && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/results ReadWrite || exit 1

# Figure out sambamba
sambamba sort -t ${THREADS} /data/bam_files/${BAM_FILE}