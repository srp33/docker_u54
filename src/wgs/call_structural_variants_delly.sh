#! /bin/bash

source usage_functions
source check_functions

set -o errexit

BAM_FILES=()
REF_GENOME=Null
OUTPUT=Null
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case "${!i}" in
    -b | --bam_file )
      check_args "${!OPTARG}" "${!i}" || exit 1
      BAM_FILES+=("/data/bam_files/${!OPTARG}")
      i=$((i+1))
      ;;
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT="${!OPTARG}"
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

[[ ${#BAM_FILES[@]} -ge 1 ]] || { echo "
ERROR: One or more BAM FILE (-b <arg>) arguments must be provided" && \
 usage_delly && exit 1; }
[[ "${REF_GENOME}" != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_delly && exit 1; }
[[ "${OUTPUT}" != "Null" ]] || { echo "
ERROR: OUTPUT (-s <arg>) argument must be provided" && \
 usage_delly && exit 1; }

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

python /check_permissions.py /data/bam_files ReadWrite || exit 1
python /check_permissions.py /data/ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"

if [[ ${REF_GENOME##*.} = "gz" ]]; then
    NEW_REF="$(echo ${REF_GENOME%.*})"
    gunzip -c /data/ref_genome/"${REF_GENOME}" > /tmp/"${NEW_REF}"
    REF_GENOME="${NEW_REF}"
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


if [[ ${VERSION_LOG} != "" ]]; then

    echo "call_somatic_variants_delly

Commands:
  delly call -g /data/ref_genome/\"${REF_GENOME}\" \"${BAM_FILES[@]}\" | \\
    bcftools view > /data/output_data/\"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python2_version )

  delly:
    version $( get_conda_version delly )

" > /data/output_data/"${VERSION_LOG}"

fi

delly call -g /tmp/"${REF_GENOME}" "${BAM_FILES[@]}" -o /tmp/"${OUTPUT%%.*}.bcf"
bcftools view /tmp/"${OUTPUT%%.*}.bcf" > /data/output_data/"${OUTPUT}"
