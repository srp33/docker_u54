#! /bin/bash

REF_GENOME=Null
THREADS=1
READ1=Null
READ2=Null

while getopts "t:r:s:h" opt; do
  case ${opt} in
    t )
      THREADS=${OPTARG}
      ;;
    r )
      REF_GENOME=${OPTARG}
      ;;
    s )
      case ${OPTARG} in
        1 )
          READ1="${!OPTIND}"
          OPTIND=$(( $OPTIND + 1 ))
          ;;
        2 )
          READ2="${!OPTIND}"
          OPTIND=$(( $OPTIND + 1 ))
          ;;
        * )
          echo "Invalid option: -s${OPTARG}"
          ;;
      esac
      ;;
    h )
      usage_align
      exit 0
      ;;
    \? )
      echo "Invalid option: -${OPTARG}" 1>&2
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
 [[ ${READ2} != "Null" ]] || { echo "
ERROR: READ 2 (-s2 <arg>) argument must be provided" && \
 usage_align && exit 1; }

# Checks for the necessary directories which are only created by volumes

[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
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

# Check for necessary index files in ref_index directory
#   If one of the files is missing, bwa index will be run

for filename in /data/ref_genome/*; do
    REF_INDEX_FILES+=($(echo "${filename##*/}"))
done

for NEEDED_FILE in ${NEEDED_FILES[@]}; do
    [[ " ${REF_INDEX_FILES[@]} " =~ " ${NEEDED_FILE} " ]] || REF_INDEXED=1
done


if [[ ${REF_INDEXED} == 1 ]]; then

    echo "The reference does not contain the proper index files. Running bwa index"

    python check_permissions.py /data/ref_genome ReadWrite || \
       { echo "Please ensure you are passing in directory and not just a file volume" && exit 1; }

    bwa index -t ${THREADS} /data/ref_genome/${REF_GENOME}
fi

bwa mem -t ${THREADS} /data/ref_genome/${REF_GENOME} \
    /data/input_data/${READ1} /data/input_data/${READ2} | \
    samtools view -@ ${THREADS} -S -b > /data/output_data/${SAMPLE}.bam

