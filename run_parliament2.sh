#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BAM_FILE=Null
REF_GENOME=Null
OUT_DIR="parliament2_out"
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -b | --bam_file )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILE="${!OPTARG}"
      i=$((i+1))
      ;;
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output )
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
      usage_delly
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
 usage_samblast && exit 1; }
[[ ${REF_GENOME} != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_align && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${REF_GENOME##*.} = "gz" ]]; then
    NEW_REF="$(echo ${REF_GENOME%.*})"
    gunzip -c /data/ref_genome/"${REF_GENOME}" > /tmp/"${NEW_REF}"
    REF_GENOME="${NEW_REF}"
    else
    ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"
fi

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && usage_parliament2 && exit 1
fi

NEEDED_INDEX=/data/ref_index/"${REF_GENOME}".fai

if [[ ! -f "${NEEDED_INDEX}" ]]; then
    echo "
    Samtools reference index (${NEEDED_INDEX}) is missing. Running samtools faidx
"
    samtools faidx /tmp/"${REF_GENOME}"
    mv /tmp/"${REF_GENOME}.fai" /data/ref_index/"${REF_GENOME}.fai"
fi

ln -s /data/ref_index/"${REF_GENOME}.fai" /tmp/"${REF_GENOME}.fai"

[[ -f /tmp/"${REF_GENOME}" ]] || { echo "
Daggummit it doesn't exist
" && exit 1; }

python /home/dnanexus/parliament2.py --bam "${BAM_FILE}" -r "${REF_GENOME}" || \
 { usage_parliament2 && exit 1; }
[[ -d /data/output_data/"${OUT_DIR}" ]] || mkdir /data/output_data/"${OUT_DIR}"
mv -r /home/dnanexus/out/* /data/output_data/"${OUT_DIR}"
