# Docker U54

## Description

This Docker image applies the `bwa mem` and `samtools view` commands to next-generation sequencing data (FASTQ files) within an isolated environment. It is designed to process data for one paired-end sample at a time. It also uses environment variables, which allow you to customize settings for these tools (for example, the reference genome that you wish to use or the number of threads you use). It assumes you have placed FASTQ files for a biological sample (two FASTQ files because they are paired-end reads) in a directory and that you wish to save the output file (BAM) for that sample to a different directory.

## Usage

When 

###### Command Reference

```
docker run -v <location of reference genome>:/data/ref_index 
-v <location of reads>:/data/sample_data
-v <location of output directory>:/data/results
--user <user id> -it [additional options]
-e REF_GENOME=<name of reference genome file>
srp33/bwamtools:latest
```

Below is an example (using example values) of how you might execute this command.

```
docker run -v /Applications/U54/ref_index:/data/ref_index 
-v /Applications/U54/in_use:/data/sample_data 
-v /Applications/U54/OutputData:/data/results 
--user 1001 
-e REF_GENOME=ucsc.hg19.fasta.gz
-it --rm srp33/bwamtools:latest
```

This command can seem rather daunting at first. Below is a description of each argument.

* `-v`
  * This flag allows the user to connect a directory on the host computer to a directory within the container. These paths are separated by a colon. Three different volumes must be specified. These are described below.
  
1. `<location of reference genome>:/data/ref_index`

   * `<location of reference genome>` should be a directory where the reference genome is
   found. This genome should have been previously indexed using `bwa index`, but if the proper
   index files are not found, the user will be prompted for either `y` or `n` (though in all
   reality, anything other than `y` will result in the container stopping) 
   as to whether the container should index the genome. Please note that this may take a 
   considerable amount of time depending on the size of the reference genome. Because of 
   possible issues with permissions (for example, if the user gives only the reference 
   genome file and not the directory, 
   resulting in the directory being hidden from the container [see `--user`]) the index files
   and a copy of the reference genome will be put in the output directory (see 
   `<location of output directory>:/data/results`). This may well be a complicated step and
   is not very efficient since it will create a copy of the reference genome, so in order to
   avoid this, __please run `bwa index` before using this container__ and __pass the entire
   directory containing the index files into the `ref_index` directory__. This will ensure that 
   index files will be found if they exist and there will be less unnecessary steps.
   
   * `/data/ref_index` is where bwamtools expects to find the reference genome. If it is not
   found, the container will fail and exit.
  
2. `<location of reads>:/data/sample_data`

   * `<location of reads>` is where the reads are located. The easiest way to use this volume
   is by putting the desired reads into a directory and simply referencing the directory
   into the volume. If there are multiple samples in the directory, the container will
   dynamically run its commands a sample at a time, outputting as many `.bam` files as there
   are samples. If the user does not want to move the files to their own directory and the user
   does not want this behavior (i.e. the user only wants to run this container on one sample)
   the user may reference the specific files as follows:
   
   `-v <location of reads>/<file name>:/data/sample_data/<file name>`
   
   * The file name does not need to coincide between the file found on the host and that on
   the container, however, the container will dynamically name the outputted `.bam` file 
   based on the name of the reads. It should also be noted that if the first part of the file
   (everything before the first `.`) does not coincide between reads of the same sample, they
   will not be analyzed as being part of the same sample. Thus consistency will allow the user
   to avoid many unforeseen issues.
   
3. `<location of output directory>:/data/results`
  
   * It is vital for the proper function of this container that the user follow two rules
   concerning the output directory:
     1. __The entire directory must be passed in__. If the user attempts to pass in a file,
     the container will not have permissions for writing within the `results` directory.
     Without these permissions, the container will fail and stop.
     2. __The docker container must have permissions for writing within the `results` directory__.
     This may done by either allowing write permissions using `chmod` although this method is not
     recommended. Instead, it is recommended that the user follow the guidelines defined in the
     `--user <user id>` section of this document.

* `--user <user id>`

  * This flag allows the user to pass read, write, and execute permissions to the container.
  Permissions are particularly important when writing, since most directories, by default,
  will block docker from writing if the proper user id is not passed in. Follow these steps
  to retrieve and use the proper user id:
  
    1. Docker does not use the user names that are probably more known to the user, but instead
    uses the user id that is generated by the computer when a user is made (without going into
    too much detail). Thus, we need to know what user id to pass in. First, let's use the 
    command `ls -lha` to view ownership of our output directory (if needed, the user 
    should create an output directory to avoid ownership conflicts):
        ```
        drwxrwxr-x  3 foo foo 4.0K Oct 16 11:21 U54
        ```
        Thus we see that `foo` owns the directory named `U54`, but we can't pass `--user foo`
        to the container since it only deals with user ids. So we must get the user id of 
        `foo`
    2. To get the user id, we use the command `id`:
        ```
        uid=1001(foo) gid=1001(foo) groups=1001(foo)
        ```
        What we are most interested in is the `uid` element (however, in this example, `foo`
        only belongs to one group, so all ids are the same, but it should be noted
        that it is likely that there will be more elements when the user runs this command,
        so just use the `uid` that coincides with the owner of the output directory).
    3. Pass the `uid` into the `docker run` command:
        ```
        docker run <volumes> --user 1001 srp33/bwamtools:latest
        ```
        This will allow the docker container to write in the output directory and avoid
        conflicts that would cause errors.
* `-it`
  * This allows the container to stop and ask for input from the user should they desire to
  run `bwa index` within the container
* `-e`
  * This allows the user to assign certain variables within the container. There are only two 
  possible options, one of which is required and the other being optional:
    1. __*REQUIRED*__ `REF_GENOME=<name of reference genome>`
       * This will allow the container to know exactly which file to use as the reference genome,
       avoiding any need for guessing or risk of getting errors. 
       * __DEFAULT:__ ucsc.hg19.fasta.gz
    2. __*OPTIONAL*__ `BWA_THREADS=<integer>`
       * The user can specify how many threads to use in `bwa mem` (and `bwa index` if it is
       run). This must be an integer.
       * __DEFAULT:__ 100
* `[additional options]`
  * Any options that the user deems necessary for `docker run`. One recommended flag is `--rm`
  since this will delete the container once it has stopped. This is useful because the
  container will automatically stop after outputting the `.bam` file and then clutter up
  Docker (not greatly, but it's still just a little messy) if not removed.
* `srp33/bwamtools:latest`
  * This should always be the last piece to the `docker run` command
  
## The Big Picture

Remember a few things:

  * You must include three volumes:
    1. `<location of reference genome>:/data/ref_index`
    2. `<location of reads>:/data/sample_data`
    3. `<location of output directory>:/data/results`
    
    Remember that each volume must be preceded by the `-v` flag and follows the format
    `<location on host>:<prescribed location in container>`. The `<location on host>` is
    determined by the user, __but the `<prescribed location in container>` must be one of the
    locations described in this document__. Otherwise the container will fail.
    
  * The docker container must have appropriate permissions:
    1. Read permissions for the reference genome directory
    2. Read permissions for the reads
    3. Read/Write permissions for the output directory
    
    Permissions are most easily dealt with using the `--user` argument, however, it is up to 
    the user to resolve any issues with permissions.
    
  * Unless you are using `ucsc.hg19.fasta.gz`, you will need to define `REF_GENOME` using the
  `-e` flag as aforementioned.
  
    