#! /bin/bash

docker run -v /Applications/U54/ref_index:/data/ref_index -v /Applications/U54/in_use:/data/SampleData zence/bwamtools:latest > output.bam
