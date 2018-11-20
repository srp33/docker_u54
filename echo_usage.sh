#! /bin/bash

echo "
Usage:
docker run \\
-v <location of reference genome>:/data/ref_genome \\
-v <location of reads>:/data/input_data \\
-v <location of output directory>:/data/output_data \\
-v <location of bam files>:/data/bam_files \\
--user <user id> -it [additional options] \\
--rm \\
srp33/somatic_wgs:latest \\
<command> <args...>

Available commands:

align
index_bam
mark_duplicates
merge_bams
slice_bam
sort_bam

See documentation for more detail
"