#! /bin/bash

#useradd -u $(id -u) dockuser

[[ $EUID -ne 0 ]] || { echo "
Cannot be run as root. Please specify user with flag --user \$(id -u):\$(id -g)
" && exit 1; }

echo "dockuser:x:$(id -u):$(id -g):,,,:/data:/bin/bash" >> /etc/passwd

#gosu dockuser -c "$@"
$@