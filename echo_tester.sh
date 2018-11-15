#! /bin/bash

source check_for_args.sh

dummy_var=()
ARGNUM=$#

for (( i=1; i<=ARGNUM; i++ )); do
  OPTARG=$((i+1))
  case ${!i} in
    -t )
      check_args "${!OPTARG}" "-t" || exit 1
      echo "${!OPTARG}"
      i=$((i+1))
      ;;
    -n )
      echo "N!!"
      ;;
    -s1 )
      echo "s1!!!"
      ;;
    -s2 )
      echo "s2!!!"
      ;;
    * )
      echo "Invalid option: ${!i}" 1>&2
      ;;
  esac
done