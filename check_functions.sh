#! /bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This command is not available to users" && exit 1
fi

check_args (){
    { [[ $1 != -* ]] && [[ $1 != "" ]]; } || { echo "$2 requires 1 argument" && exit 1; }
}

get_conda_version (){
    conda list --json $1 | python -c "import sys, json; print(json.load(sys.stdin)[0]['version'])"
}

get_python_version (){
    python -c "import sys; print('.'.join(str(x) for x in sys.version_info[:3]))"
}

get_python2_version (){
    python2.7 -c "import sys; print('.'.join(str(x) for x in sys.version_info[:3]))"
}