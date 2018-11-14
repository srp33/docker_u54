#! /bin/bash

BAM_FILE=Null
REGION=Null
OUTPUT=Null
THREADS=1

while getopts "b:t:r:o:h" opt; do
  case ${opt} in
    b )
      BAM_FILE=${OPTARG}
      ;;
    t )
      THREADS=${OPTARG}
      ;;
    r )
      REGION=${OPTARG}
      ;;
    o )
      OUTPUT=${OPTARG}
      ;;
    h )
      usage_slice_bam
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
 usage_slice_bam && exit 1; }
[[ ${REGION} != "Null" ]] || { echo "
ERROR: REGION (-r <arg>) argument must be provided" && \
 usage_slice_bam && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_slice_bam && exit 1; }

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

sambamba slice -o /data/output_data/${OUTPUT} /data/bam_files/${BAM_FILE} ${REGION}