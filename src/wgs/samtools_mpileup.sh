#! /bin/bash

source usage_functions
source check_functions

REF_GENOME=Null
NORMAL=Null
TUMOR=Null
OUTPUT=Null
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME=${!OPTARG}
      i=$((i+1))
      ;;
    -n | --normal )
      check_args "${!OPTARG}" "${!i}" || exit 1
      NORMAL="${!OPTARG}"
      i=$((i+1))
      ;;
    -t | --tumor )
      check_args "${!OPTARG}" "${!i}" || exit 1
      TUMOR="${!OPTARG}"
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
      usage_mpileup
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done


REF_INDEX_FILES=()
NEEDED_FILES=(${REF_GENOME}.amb ${REF_GENOME}.ann ${REF_GENOME}.bwt ${REF_GENOME}.pac ${REF_GENOME}.sa)
REF_INDEXED=0
MISSING_VOLUMES=()
EXIT_CODE=0

# Since the default value for REF_GENOME and SAMPLE is "Null", we can check if those have been changed

[[ ${REF_GENOME} != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_mpileup && exit 1; }
[[ ${NORMAL} != "Null" ]] || { echo "
ERROR: NORMAL (-n <arg>) argument must be provided" && \
 usage_mpileup && exit 1; }
[[ ${TUMOR} != "Null" ]] || { echo "
ERROR: TUMOR (-t <arg>) argument must be provided" && \
 usage_mpileup && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_mpileup && exit 1; }

# Checks for the necessary directories which are only created by volumes

[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

# Check permissions of each directory

python /check_permissions.py /data/ref_genome Read ${REF_GENOME} || exit 1
python /check_permissions.py /data/bam_files Read ${NORMAL} || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1

#mkdir /tmp
ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"

# Check for necessary index files in ref_index directory
#   If one of the files is missing, bwa index will be run

INDEX=$(echo ${REF_GENOME} | grep -o '\.' | grep -c '\.')
if [[ ${REF_GENOME: -${INDEX}} = ".gz" ]]; then
    NEW_REF="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})"
    gunzip -c /data/ref_genome/"${REF_GENOME}" > /tmp/"${NEW_REF}"
    REF_GENOME="${NEW_REF}"
fi

NEEDED_FILE=/data/ref_index/${REF_GENOME}.fai

if [[ ! -f "${NEEDED_FILE}" ]]; then
    echo "
    Samtools reference index (${NEEDED_FILE}) is missing. Running samtools faidx
"
    samtools faidx /tmp/"${REF_GENOME}"
    mv /tmp/"${REF_GENOME}.fai" /data/ref_index/"${REF_GENOME}.fai"
fi

ln -s /data/ref_index/"${REF_GENOME}.fai" /tmp/"${REF_GENOME}.fai"

if [[ ${VERSION_LOG} != "" ]]; then

    echo "bwa_mem_align

Command:
  samtools mpileup -B -q 1 -f /tmp/"${REF_GENOME}" \\
    /data/bam_files/\"${NORMAL}\" /data/bam_files/\"${TUMOR}\" > \"${OUTPUT}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  samtools:
    version $( get_conda_version samtools )
" > /data/output_data/"${VERSION_LOG}"

fi

samtools mpileup -B -q 1 -f /tmp/"${REF_GENOME}" \
/data/bam_files/"${NORMAL}" /data/bam_files/"${TUMOR}" > /data/output_data/"${OUTPUT}"