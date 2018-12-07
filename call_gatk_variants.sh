#! /bin/bash

source usage_functions
source check_functions

TUMOR=Null
NORMAL=Null
OUTPUT=Null
REF_GENOME=Null
TUMOR_SAMPLE=Null
NORMAL_SAMPLE=Null
VERSION_LOG=""
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
    -ts | --tumor-sample )
      check_args "${!OPTARG}" "${!i}" || exit 1
      TUMOR_SAMPLE=${!OPTARG}
      i=$((i+1))
      ;;
    -ns | --normal-sample )
      check_args "${!OPTARG}" "${!i}" || exit 1
      NORMAL_SAMPLE=${!OPTARG}
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
    --log )
      check_args "${!OPTARG}" "${!i}" || exit 1
      VERSION_LOG=${!OPTARG}
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
ERROR: TUMOR BAM FILE (-t <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${TUMOR_SAMPLE} != "Null" ]] || { echo "
ERROR: TUMOR SAMPLE (-ts <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${NORMAL} != "Null" ]] || { echo "
ERROR: NORMAL BAM FILE (-n <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${NORMAL_SAMPLE} != "Null" ]] || { echo "
ERROR: NORMAL SAMPLE (-n <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${OUTPUT} != "Null" ]] || { echo "
ERROR: OUTPUT (-o <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${REF_GENOME} != "Null" ]] || { echo "
ERROR: REFERENCE GENOME (-r <arg>) argument must be provided" && \
 usage_mutect && exit 1; }
[[ ${TUMOR_SAMPLE} != "Null" ]] || TUMOR_SAMPLE="${TUMOR%%.*}"
[[ ${NORMAL_SAMPLE} != "Null" ]] || NORMAL_SAMPLE="${NORMAL%%.*}"

EXIT_CODE=0
REF_LOCATION=/data/ref_genome/${REF_GENOME}
NEEDED_INDEX=/data/ref_genome/${REF_GENOME}.fai
MISSING_VOLUMES=()

[[ -d /data/bam_files ]] || { MISSING_VOLUMES+=(/data/bam_files) && EXIT_CODE=1; }
[[ -d /data/ref_genome ]] || { MISSING_VOLUMES+=(/data/ref_genome) && EXIT_CODE=1; }
[[ -d /data/output_data ]] || { MISSING_VOLUMES+=(/data/output_data) && EXIT_CODE=1; }

if [[ ${EXIT_CODE} = 1 ]]; then
    echo "
    The following volumes are missing: ${MISSING_VOLUMES[@]}" && echo_usage && exit 1
fi

if [[ ! -f ${NEEDED_INDEX} ]]; then
    echo "
    Samtools reference index (${NEEDED_INDEX}) is missing. Running samtools faidx
"
    if [[ ${REF_GENOME: -3} = ".gz" ]]; then
        mkdir /data/temp_ref
        INDEX=$(echo ${REF_GENOME} | grep -o '\.' | grep -c '\.')
        NEW_REF="$(echo ${REF_GENOME} | cut -d '.' -f -${INDEX})"
        gunzip -c /data/ref_genome/${REF_GENOME} > /data/temp_ref/${NEW_REF}
        REF_LOCATION=/data/temp_ref/${NEW_REF}
    fi
    samtools faidx ${REF_LOCATION}
fi

INDEX=$(echo ${REF_LOCATION} | grep -o '\.' | grep -c '\.')
NEEDED_DICT="$(echo ${REF_LOCATION} | cut -d '.' -f -${INDEX})".dict

if [[ ! -f ${NEEDED_DICT} ]]; then
    echo "
    Fasta dict file (${NEEDED_DICT}) is missing. Running gatk CreateSequenceDictionary
"
    gatk CreateSequenceDictionary --REFERENCE "${REF_LOCATION}"
fi

python /check_permissions.py /data/bam_files Read "${TUMOR}" || exit 1
python /check_permissions.py /data/ref_genome Read "${REF_GENOME}" || exit 1
python /check_permissions.py /data/output_data ReadWrite || exit 1

if [[ ${VERSION_LOG} != "" ]]; then

    echo "call_gatk_variants

Commands:
  gatk Mutect2 -I /data/bam_files/\"${TUMOR}\" -tumor \"${TUMOR_SAMPLE}\" \\
    -I /data/bam_files/"${NORMAL}" -normal \"${NORMAL_SAMPLE}\" \\
    -O /data/results/\"${OUTPUT}\" -R \"${REF_LOCATION}\"

Timestamp: $(date '+%d/%m/%Y %H:%M:%S')

Software used:
  Bash:
    $( bash --version )

  Python:
    version $( get_python_version )

  gatk:
    version 4.0.11.0
" > /data/output_data/"${VERSION_LOG}"

fi

gatk Mutect2 -I /data/bam_files/"${TUMOR}" -tumor "${TUMOR_SAMPLE}" \
  -I /data/bam_files/"${NORMAL}" -normal "${NORMAL_SAMPLE}" \
  -O /data/results/"${OUTPUT}" -R "${REF_LOCATION}"