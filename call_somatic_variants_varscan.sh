#! /bin/bash

source usage_functions
source check_functions

PILEUP=Null
OUTPUT=Null
REF_GENOME=""
NORMAL=""
TUMOR=""
RUN_PILEUP=0
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -p | --pileup )
      check_args "${!OPTARG}" "${!i}" || exit 1
      PILEUP=${!OPTARG}
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT=${!OPTARG}
      i=$((i+1))
      ;;
    --run_pileup )
      RUN_PILEUP=1
      ;;
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME=${!OPTARG}
      i=$((i+1))
      ;;
    -n | --normal )
      check_args "${!OPTARG}" "${!i}" || exit 1
      NORMAL=${!OPTARG}
      i=$((i+1))
      ;;
    -t | --tumor )
      check_args "${!OPTARG}" "${!i}" || exit 1
      TUMOR=${!OPTARG}
      i=$((i+1))
      ;;
    -h | --help )
      usage_strelka
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${PILEUP} != "Null" ]] || { echo "
ERROR: PILEUP FILE (-p <arg>) argument must be provided" && \
 usage_varscan && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_varscan && exit 1; }


EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/input_data ]] || { MISSING_VOLUMES+=(/data/input_data) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${RUN_PILEUP} -eq 1 ]]; then
    echo "
Running samtools_mpileup
"
    samtools_mpileup -t "${TUMOR}" -n "${NORMAL}" -r "${REF_GENOME}" -o "${PILEUP}"
    ln -s /data/output_data/"${PILEUP}" /tmp/"${PILEUP}"
else
    python /check_permissions.py /data/input_data Read "${PILEUP}" || exit 1
    ln -s /data/input_data/"${PILEUP}" /tmp/"${PILEUP}"
fi

varscan somatic /tmp/"${PILEUP}" /data/output_data/"${OUTPUT}" --mpileup 1