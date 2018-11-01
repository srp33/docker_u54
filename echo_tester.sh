#! /bin/bash

while getopts "t:n:" opt; do
  case ${opt} in
    t )
      threads="$OPTARG"
      ;;
    n )
      dummy_var="$OPTARG"
      ;;
    \? )
      echo "Invalid option: -$OPTARG" 1>&2
      ;;
  esac
done

echo ${threads}
echo ${dummy_var}