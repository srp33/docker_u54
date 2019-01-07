#! /bin/bash

echo "Additional details about this software can be found at https://github.com/srp33/somatic_wgs

Available commands:
"
ls /usr/local/bin/wgs | awk "{print $1;}"