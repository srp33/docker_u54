#! /bin/bash

source usage_functions
source check_for_args

TUMOR=Null
NORMAL=Null
OUTPUT=Null
REF_GENOME=Null
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -t | --tumor )
      check_args "${!OPTARG}" "${!i}" || exit 1
      TUMOR=${!OPTARG}
      i=$((i+1))
      ;;
    -n | --normal )
      check_args "${!OPTARG}" "${!i}" || exit 1
      NORMAL=${!OPTARG}
      i=$((i+1))
      ;;
    -o | --output )
      check_args "${!OPTARG}" "${!i}" || exit 1
      OUTPUT=${!OPTARG}
      i=$((i+1))
      ;;
    -r | --reference )
      check_args "${!OPTARG}" "${!i}" || exit 1
      REF_GENOME=${!OPTARG}
      i=$((i+1))
      ;;
    -h | --help )
      usage_mutect
      exit 0
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      exit 1
      ;;
  esac
done

[[ ${TUMOR} != "Null" ]] || { echo "
ERROR: TUMOR BAM FILE (-p <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${NORMAL} != "Null" ]] || { echo "
ERROR: NORMAL BAM FILE (-t <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_mutect && exit 1; }

gatk MuTect2 -I:tumor /data/bam_files/${TUMOR} -I:normal /data/bam_files/${NORMAL} \
  -O /data/results/${OUTPUT} -R /data/ref_genome/${REF_GENOME}