#! /bin/bash

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