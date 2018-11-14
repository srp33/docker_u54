#! /bin/bash

dummy_var=()

while getopts "t:n:s:" opt; do
  case ${opt} in
    t )
      threads="$OPTARG"
      ;;
    n )
      dummy_var+=("/bin/bash/${OPTARG}")
      ;;
    s )
      case "${OPTARG}" in
        1 )
          val="${!OPTIND}"
          echo ${val}
          echo ${OPTIND}
          echo "woah buddy ${OPTARG}"
          ;;
        2 )
          echo "noah buddy"
          ;;
        * )
          echo "Invalid option: -s${OPTARG}"
          ;;
      esac;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      ;;
  esac
done