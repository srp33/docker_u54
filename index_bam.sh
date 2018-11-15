#! /bin/bash

source check_for_args

BAM_FILE=Null
THREADS=1
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -b )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE=${OPTARG}
      i=$((i+1))
      ;;
    -t )
      check_args "${!OPTARG}" "${!i}" || exit 1
      THREADS=${OPTARG}
      i=$((i+1))
      ;;
    -h )
      usage_index_bam
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
 usage_index_bam && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1


sambamba index -t ${THREADS} /data/bam_files/${BAM_FILE} /data/bam_files/${BAM_FILE}.bai