#! /bin/bash

source usage_functions
source check_functions

set -o errexit

REF_GENOME=Null
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG="${!OPTARG}"
      i=$((i+1))
      ;;
    -h | --help )
      usage_index_fasta
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${REF_GENOME} != "Null" ]] || { echo "
ERROR: REF_GENOME (-r <arg>) argument must be provided" && \
 usage_index_fasta && exit 1; }

EXIT_CODE=0
MISSING_VOLUMES=()

[[ -d ref_genome ]] || { MISSING_VOLUMES+=(ref_genome) && EXIT_CODE=1; }
[[ -d output_data ]] || { MISSING_VOLUMES+=(output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

# Check permissions of each directory

python /check_permissions.py ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

  echo "index_fasta

Commands:
  bwa index -a bwtsw /tmp/"${REF_GENOME}"
  samtools faidx /tmp/"${REF_GENOME}"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  bwa:
    version $( get_conda_version bwa )

  samtools:
    version $( get_conda_version samtools )
" > /data/output_data/"${VERSION_LOG}"
fi

ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"

bwa index -a bwtsw /tmp/"${REF_GENOME}"
samtools faidx /tmp/"${REF_GENOME}"

REF_INDEX_FILES=("${REF_GENOME}".amb "${REF_GENOME}".ann "${REF_GENOME}".bwt "${REF_GENOME}".pac "${REF_GENOME}".sa "${REF_GENOME}".fai)

for REF_INDEX_FILE in ${REF_INDEX_FILES[@]}; do
    mv /tmp/"${REF_INDEX_FILE}" /data/output_data
done
