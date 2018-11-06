#! /bin/bash

echo "Usage:
docker run \\
-v <location of reference genome>:/data/ref_index \\
-v <location of reads>:/data/sample_data \\
-v <location of output directory>:/data/results \\
-v <location of bam files>:/data/bam_files \\
--user <user id> -it [additional options] \\
--rm \\
srp33/u54:latest \\
<command> <args...>"