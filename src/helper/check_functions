#! /bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This command is not available to users" && exit 1
fi

check_args (){
    { [[ $1 != -* ]] && [[ $1 != "" ]]; } || { echo "$2 requires 1 argument" && exit 1; }
}

get_bash_version (){
    bash --version
    echo
}

get_conda_version (){
    echo "${1}:"
    conda list --json $1 | python -c "import sys, json; print(json.load(sys.stdin)[0]['version'])"
    echo
}

get_python_version (){
    echo "Python:"
    python -c "import sys; print('.'.join(str(x) for x in sys.version_info[:3]))"
    echo
}

get_python2_version (){
    echo "Python:"
    python2.7 -c "import sys; print('.'.join(str(x) for x in sys.version_info[:3]))"
    echo
}

get_file_extension (){
    filename=$(basename -- "$1")
    extension="${filename##*.}"
    echo "$extension"
}

get_file_without_extension (){
    extension="$(get_file_extension $1)"
    basename "$1" ".${extension}"
}
