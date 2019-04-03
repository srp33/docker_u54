#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BAM_FILE=Null
REF_GENOME=Null
THREADS=1
OUT_DIR=Null
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -b | --bam_file )
      check_args "${!OPTARG}" "{$i}" || exit 1
      BAM_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -r | --reference_genome )
      check_args "${!OPTARG}" "{$i}" || exit 1
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output_directory )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUT_DIR="${!OPTARG}"
      i=$((i+1))
      ;;
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_sve
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ "${BAM_FILE}" != "Null" ]] || { echo "
ERROR: BAM FILE (-b <arg>) argument must be provided" && \
 usage_sve && exit 1; }
[[ "${REF_GENOME}" != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_sve && exit 1; }
[[ "${OUT_DIR}" != "Null" ]] || { echo "
ERROR: OUTPUT DIRECTORY (-o <arg>) argument must be provided" && \
 usage_sve && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

python /check_permissions.py /data/bam_files Read ${BAM_FILE} || exit 1
python /check_permissions.py /data/ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

    echo "run_sve

Commands:
  /usr/local/bin/python2.7 /tools/SVE/scripts/FusorSV/FusorSV.py \
  -r /tmp/\"${REF_GENOME}\" \
  -i /data/input_data/\"${BAM_FILE}\" \
  -p \"${THREADS}\" \
  -o /data/output_data/\"${OUT_DIR}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( /usr/local/bin/python2.7 -c "import sys; print('.'.join(str(x) for x in sys.version_info[:3]))" )

  FusorSV:
    version 0.1.0

" > /data/output_data/"${VERSION_LOG}"

fi

SV_CALLERS=("breakdancer" "breakseq" "cnvnator" "hydra" "delly" "lumpy" "cnmops")

for caller in ${SV_CALLERS[@]}; do
  timeout 6h /usr/local/bin/python2.7 /tools/SVE/bin/sve call \
                -r /data/ref_genome/${REF_GENOME} \
                -o /data/output_data/${OUT_DIR} \
                -a ${caller} \
                -g hg19 \
                /data/bam_files/${BAM_FILE}
done