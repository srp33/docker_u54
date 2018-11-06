#! /bin/bash

BAM_FILE=Null
REGION=Null
OUTPUT=Null
THREADS=1

while getopts "b:t:r:o:" opt; do
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
    \? )
      echo "Invalid option: -${OPTARG}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${BAM_FILE} != "Null" ]] || { echo "ERROR: BAM FILE (-b <arg>) argument must be provided" && exit 1; }
[[ ${REGION} != "Null" ]] || { echo "ERROR: REGION (-r <arg>) argument must be provided" && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "ERROR: OUTPUT (-o <arg>) argument must be provided" && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/results ]] || { MISSING_VOLUMES+=(/data/results) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/results ReadWrite || exit 1

sambamba slice -o /data/results/${OUTPUT} /data/bam_files/${BAM_FILE} ${REGION}