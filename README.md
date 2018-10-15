# Docker U54

## Description

* This container performs the `bwa mem` and `samtools view` commands within an isolated environment

## Usage

The command for this container can seem rather daunting at first:

`docker run -v /Applications/U54/ref_index:/data/ref_index 
-v /Applications/U54/in_use:/data/sample_data 
-v /Applications/U54/OutputData:/data/results 
--user 1004 srp33/bwamtools:latest`

so I will step through each piece of this command so using it will be smooth and simple.

* `-v`
  * This flag allows you to connect a directory or folder to an identical directory or folder
  in the container. There are three important volumes that are required for proper use of the
  bwamtools container
  
1. `<location of reference genome>:/data/ref_index`

   * `<location of reference genome>` should be a directory where the reference genome is
   found. This genome should have been previously indexed using `bwa index`, but if the proper
   index files are not found, the user will be prompted for either `y` or `n` (though in all
   reality, anything other than `y` will result in a failure) as to whether the container
   should index the genome. Please note that this may take a considerable amount of time 
   depending on the size of the reference genome.
   
   * `/data/ref_index` is where bwamtools expects to find the reference genome. If it is not
   found, the container will fail and exit.
  
2. `<location of reads>:/data/sample_data`

   * `<location of reads>` is where the reads are located. The easiest way to use this volume
   is by putting the desired reads into a directory and simply referencing the directory
   into the volume. If there are multiple samples in the directory, the container will
   dynamically run its commands a sample at a time, outputting as many `.bam` files as there
   are samples. If you do not want to move the files to their own directory and you do not
   want this behavior (i.e. you only want to run this container on one sample) you may 
   reference the specific files as follows:
   
   `-v <location of reads>/<file name>:/data/sample_data/<file name>`
   
   * The file name does not need to coincide between the file found on the host and that on
   the container, however, the container will dynamically name the outputted `.bam` file 
   based on the name of the reads. It should also be noted that if the first part of the file
   (everything before the first `.`) does not coincide between reads of the same sample, they
   will not be analyzed as being part of the same sample. Thus consistency will allow the user
   to avoid many unforeseen issues.