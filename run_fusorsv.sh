#! /bin/bash

source usage_functions
source check_functions

set -o errexit

VCF_FILES=()
REF_GENOME=""

/usr/local/bin/python2.7 /tools/SVE/scripts/FusorSV/FusorSV.py -h