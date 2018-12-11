#! /bin/bash

source usage_functions
source check_functions

REF_GENOME=Null
THREADS=1
READ1=Null
READ2=""
OUTPUT=Null
VERSION_LOG=""
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -t | --nthreads )
      check_args "${!OPTARG}" "${!i}" || exit 1
      THREADS=${!OPTARG}
      i=$((i+1))
      ;;
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME=${!OPTARG}
      i=$((i+1))
      ;;
    -s1 | --sample1 )
      check_args "${!OPTARG}" "${!i}" || exit 1
      READ1="${!OPTARG}"
      i=$((i+1))
      ;;
    -s2 | --sample2 )
      check_args "${!OPTARG}" "${!i}" || exit 1
      READ2="<(zcat /data/input_data/\"${!OPTARG}\" | awk 'int((NR-1)/4)%7==6')"
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
      usage_align
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
 usage_align && exit 1; }
[[ ${READ1} != "Null" ]] || { echo "
ERROR: READ 1 (-s1 <arg>) argument must be provided" && \
 usage_align && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_align && exit 1; }

# Checks for the necessary directories which are only created by volumes

[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/ref_index ]] || { MISSING_VOLUMES+=(/data/ref_index) && EXIT_CODE=1; }
[[ -d /data/input_data ]] || { MISSING_VOLUMES+=(/data/input_data) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

# Check permissions of each directory

python /check_permissions.py /data/ref_genome Read ${REF_GENOME} || exit 1
python /check_permissions.py /data/input_data Read ${READ1} || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1
python /check_permissions.py /data/ref_index ReadWrite || exit 1

#mkdir /tmp
ln -s /data/ref_genome/"${REF_GENOME}" /tmp/"${REF_GENOME}"

# Check for necessary index files in ref_index directory
#   If one of the files is missing, bwa index will be run

for filename in /data/ref_index/*; do
    REF_INDEX_FILES+=($(echo "${filename##*/}"))
done

for NEEDED_FILE in ${NEEDED_FILES[@]}; do
    { [[ " ${REF_INDEX_FILES[@]} " =~ " ${NEEDED_FILE} " ]] \
&& ln -s /data/ref_index/"${NEEDED_FILE}" /tmp/"${NEEDED_FILE}"; } || REF_INDEXED=1
done


if [[ ${REF_INDEXED} == 1 ]]; then

    echo "The reference does not contain the proper index files. Running bwa index"

    bwa index /tmp/"${REF_GENOME}"

    for NEEDED_FILE in ${NEEDED_FILES[@]}; do
        mv /tmp/"${NEEDED_FILE}" /data/ref_index
        ln -s /data/ref_index/"${NEEDED_FILE}" /tmp/"${NEEDED_FILE}"
    done
fi

if [[ ${VERSION_LOG} != "" ]]; then

    if [[ ${READ2} = "" ]]; then
        COMMAND="bwa mem -t ${THREADS} /tmp/\"${REF_GENOME}\" \\
    <(zcat /data/input_data/\"${READ1}\" | awk 'int((NR-1)/4)%7==6') |\\
    samtools view -@ ${THREADS} -S -b > /data/output_data/\"${OUTPUT}\""
    else
        COMMAND="bwa mem -t ${THREADS} /tmp/\"${REF_GENOME}\" \\
    <(zcat /data/input_data/\"${READ1}\" | awk 'int((NR-1)/4)%7==6') \\
    ${READ2} |\\
    samtools view -@ ${THREADS} -S -b > /data/output_data/\"${OUTPUT}\""
    fi

    echo "bwa_mem_align

Command:
  ${COMMAND}

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

if [[ ${READ2} = "" ]]; then
    bwa mem -t ${THREADS} /tmp/"${REF_GENOME}" \
        <(zcat /data/input_data/"${READ1}" | awk 'int((NR-1)/4)%7==6') |\
        samtools view -@ ${THREADS} -S -b > /data/output_data/"${OUTPUT}"
else
    bwa mem -t ${THREADS} /tmp/"${REF_GENOME}" \
        <(zcat /data/input_data/"${READ1}" | awk 'int((NR-1)/4)%7==6') \
        "${READ2}" |\
        samtools view -@ ${THREADS} -S -b > /data/output_data/"${OUTPUT}"

fi
