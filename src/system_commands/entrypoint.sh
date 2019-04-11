#! /bin/bash

[[ $EUID -ne 0 ]] || { echo "
Cannot be run as root. Please specify user with flag --user \$(id -u):\$(id -g)
" && exit 1; }

$@