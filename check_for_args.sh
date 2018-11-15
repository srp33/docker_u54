#! /bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This command is not available to users" && exit 1
fi

check_args (){
    { [[ $1 != -* ]] && [[ $1 != "" ]]; } || { echo "$2 requires 1 argument" && exit 1; }
}