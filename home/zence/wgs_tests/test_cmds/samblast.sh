#! /bin/bash

docker run \\
  -v /Applications/U54/bam_files:/data/bam_files \\
  -v /Applications/:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  samblast \\
    -b 101024.zence.bam