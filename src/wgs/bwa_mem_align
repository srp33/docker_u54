#! /bin/bash

source usage_functions
source check_functions

set -o errexit

REF_GENOME=Null
THREADS=1
READ1=Null
READ2=""
OUTPUT=Null
VERSION_LOG=""
CHUNKS=1
PROCESS_CHUNK=0
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
      REF_GENOME="${!OPTARG}"
      i=$((i+1))
      ;;
    -s1 | --sample1 )
      check_args "${!OPTARG}" "${!i}" || exit 1
      READ1="${!OPTARG}"
      i=$((i+1))
      ;;
    -s2 | --sample2 )
      check_args "${!OPTARG}" "${!i}" || exit 1
      READ2="${!OPTARG}"
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT="${!OPTARG}"
      i=$((i+1))
      ;;
    -c | --nchunks )
      check_args "${!OPTARG}" "${!i}" || exit 1
      CHUNKS="${!OPTARG}"
      i=$((i+1))
      ;;
    -p | --process_chunk )
      check_args "${!OPTARG}" "${!i}" || exit 1
      PROCESS_CHUNK="${!OPTARG}"
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
NEEDED_FILES=("${REF_GENOME}".amb "${REF_GENOME}".ann "${REF_GENOME}".bwt "${REF_GENOME}".pac "${REF_GENOME}".sa)
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
[[ ${CHUNKS} -ge 1 ]] || { echo "
ERROR: CHUNKS [${CHUNKS}](-c <arg>) argument must be greater than or equal to 1" && \
 usage_align && exit 1; }
[[ ${PROCESS_CHUNK} -lt ${CHUNKS} ]] || { echo "
ERROR: PROCESS_CHUNK [${PROCESS_CHUNK}] (-p <arg>) argument must be smaller than \
CHUNKS [${CHUNKS}] (-c <arg>) argument" && \
 usage_align && exit 1; }

# Checks for the necessary directories which are only created by volumes

[[ -d ref_genome ]] || { MISSING_VOLUMES+=(ref_genome) && EXIT_CODE=1; }
[[ -d ref_index ]] || { MISSING_VOLUMES+=(ref_index) && EXIT_CODE=1; }
[[ -d input_data ]] || { MISSING_VOLUMES+=(input_data) && EXIT_CODE=1; }
[[ -d output_data ]] || { MISSING_VOLUMES+=(output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

# Check permissions of each directory

python /check_permissions.py ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py input_data Read "${READ1}" || exit 1
python /check_permissions.py output_data ReadWrite || exit 1
python /check_permissions.py ref_index ReadWrite || exit 1

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

if [[ ${CHUNKS} -gt 1 ]]; then
    INDEX=$(echo ${OUTPUT} | grep -o '\.' | grep -c '\.')
    OUT_NAME="$(echo ${OUTPUT} | cut -d '.' -f -${INDEX})"
fi

if [[ ${REF_INDEXED} == 1 ]]; then

    echo "The reference does not contain the proper index files. Running bwa index"

    bwa index /tmp/"${REF_GENOME}"

    for NEEDED_FILE in ${NEEDED_FILES[@]}; do
        mv /tmp/"${NEEDED_FILE}" ref_index
        ln -s /data/ref_index/"${NEEDED_FILE}" /tmp/"${NEEDED_FILE}"
    done
fi

if [[ ${VERSION_LOG} != "" ]]; then

    if [[ ${READ2} = "" ]]; then
        if [[ ${CHUNKS} -gt 1 ]]; then
            COMMAND="bwa mem -t ${THREADS} /tmp/\"${REF_GENOME}\" \\
    <(zcat /data/input_data/\"${READ1}\" | awk -v chunks=${CHUNKS} -v process_chunk=${PROCESS_CHUNK} \\
      'int((NR-1)/4)%chunks==process_chunk') |\\
    samtools view -@ ${THREADS} -S -b > /data/output_data/\"${OUT_NAME}\".${PROCESS_CHUNK}.bam"
        else
            COMMAND="bwa mem -t ${THREADS} /tmp/\"${REF_GENOME}\" \
    /data/input_data/\"${READ1}\" | \
    samtools view -@ ${THREADS} -S -b > /data/output_data/\"${OUTPUT}\""
        fi
    else
        if [[ ${CHUNKS} -gt 1 ]]; then
            COMMAND="bwa mem -t ${THREADS} /tmp/\"${REF_GENOME}\" \\
    <(zcat /data/input_data/\"${READ1}\" | awk -v chunks=${CHUNKS} -v process_chunk=${PROCESS_CHUNK} \\
      'int((NR-1)/4)%chunks==process_chunk') \\
    <(zcat /data/input_data/\"${READ2}\" | awk -v chunks=${CHUNKS} -v process_chunk=${PROCESS_CHUNK} \\
      'int((NR-1)/4)%chunks==process_chunk') | \\
    samtools view -@ ${THREADS} -S -b > /data/output_data/\"${OUT_NAME}\".${PROCESS_CHUNK}.bam"
        else
            COMMAND="bwa mem -t ${THREADS} /tmp/\"${REF_GENOME}\" \\
    /data/input_data/\"${READ1}\" /data/input_data/\"${READ2}\" | \\
    samtools view -@ ${THREADS} -S -b > /data/output_data/\"${OUTPUT}\""
        fi
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
    if [[ ${CHUNKS} -gt 1 ]]; then
        bwa mem -t ${THREADS} /tmp/"${REF_GENOME}" \
            <(zcat /data/input_data/"${READ1}" | awk -v chunks=${CHUNKS} -v process_chunk=${PROCESS_CHUNK} \
              'int((NR-1)/4)%chunks==process_chunk') | \
            samtools view -@ ${THREADS} -S -b > /data/output_data/"${OUT_NAME}".${PROCESS_CHUNK}.bam
    else
        #if [[ ${WORK_QUEUE} = 1 ]]; then
        #    use_make_queue \
        #        -c "bwa" \
        #        -c "samtools" \
        #        -i /data/input_data/"${READ1}" \
        bwa mem -t ${THREADS} /tmp/"${REF_GENOME}" \
            /data/input_data/"${READ1}" | \
            samtools view -@ ${THREADS} -S -b > /data/output_data/"${OUTPUT}"
    fi
else
    if [[ ${CHUNKS} -gt 1 ]]; then
        bwa mem -t ${THREADS} /tmp/"${REF_GENOME}" \
            <(zcat /data/input_data/"${READ1}" | awk -v chunks=${CHUNKS} -v process_chunk=${PROCESS_CHUNK} \
              'int((NR-1)/4)%chunks==process_chunk') \
            <(zcat /data/input_data/"${READ2}" | awk -v chunks=${CHUNKS} -v process_chunk=${PROCESS_CHUNK} \
              'int((NR-1)/4)%chunks==process_chunk') | \
            samtools view -@ ${THREADS} -S -b > /data/output_data/"${OUT_NAME}".${PROCESS_CHUNK}.bam
    else
        bwa mem -t ${THREADS} /tmp/"${REF_GENOME}" \
            /data/input_data/"${READ1}" /data/input_data/"${READ2}" | \
            samtools view -@ ${THREADS} -S -b > /data/output_data/"${OUTPUT}"
    fi

fi
