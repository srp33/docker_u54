# Framework for performing whole-genome variant calling across dispersed computing environments

## Description

This repository contains scripts that enable researchers to manipulate BAM files and perform whole-genome variant 
calling on Illumina sequencing data. We package software within a Docker image and enable the 
user to execute the software via commands and arguments. We have also configured the image so that data can easily be passed into and out of a Docker container.

For a high-level overview of how Docker works (and software containers in general), read [this paper](https://gigascience.biomedcentral.com/articles/10.1186/s13742-016-0135-4).

To get started, install the [Docker CE](https://docs.docker.com/install) software on your computer.

Currently, we support the following commands:

* `bwa_mem_align` (align FASTA files to a reference genome using `bwa mem`)
* `sort_bam` (sort a BAM file using `sambamba`)
* `index_bam` (index a BAM file using `sambamba`)
* `mark_duplicates` (mark duplicates in a BAM file using `sambamba`)
* `slice_bam` (slice/split a BAM file using `sambamba`)
* `merge_bams` (merge BAM files using `sambamba`)
* `add_read_groups` (add read groups to a BAM file)
* `samtools_mpileup.sh` (create pileup file using `samtools`)
* `call_somatic_variants_gatk4` (call variants using Mutect2 (GATK 4))
* `call_somatic_variants_strelka` (call variants using Strelka)
* `call_structural_variants_manta` (call variants using Manta)
* `call_structural_variants_lumpy` (call variants using Lumpy)
* `call_structural_variants_delly` (call variants using Delly2)
* `base_recalibrator` (detect systematic errors in base quality scores (GATK 4))
* `apply_bqsr` (apply base quality score recalibration (GATK 4))
* `samblast` (create split and discordant read SAM files from BAM file using `samblaster`)
* `run_survivor` (merge VCF files created through structural variant calling (SURVIVOR))
* `svtype_vcf` (genotype a VCF file (SVTyper))

The default behavior of the Docker container is to display the usage and available commands. If the user executes the following command, usage information will be displayed.

```
docker run --rm srp33/somatic_wgs:latest
```

For any of the commands listed above, the user can specify the `-h` flag to view available options for 
that command. For example, the following command would provide usage information for the 
`bwa_mem_align` command, which uses [bwa mem](https://github.com/lh3/bwa) to align sequencing reads to a reference genome.

```
docker run --rm srp33/somatic_wgs:latest bwa_mem_align -h
```

Below is an example of how the `bwa_mem_align` command could be executed.

```
docker run \
  -v /MyData/Reference_Genomes/hg19:/data/ref_genome \ 
  -v /MyData/FASTQ:/data/input_data \
  -v /MyData/BAM:/data/output_data \
  --user $(id -u):$(id -g) \
  --rm \
  srp33/somatic_wgs:latest \
  bwa_mem_align \
    -r ucsc.hg19.fasta.gz \
    -s1 101024.1.fastq.gz \
    -s2 101024.2.fastq.gz \
    -o 101024.bam \
    -t 10
```

A variety of arguments must be specified, as described below:

* The first three arguments (each beginning with `-v`) specify [volumes](https://docs.docker.com/storage/volumes). Volumes enable data to be shared between the host operating system and the Docker container. The path specified before each colon indicates a directory on the host; the path specified after the colon indicates the corresponding directory within the container (these are static; please do not change them).
    - In the example above, the first volume specifies the location of the reference genome. This should be a directory that contains a FASTA file (can be gzipped) and the index for the reference genome. If the index does not already exist, our scripts will create it using `samtools`.
    - The second volume specifies the directory where input files are stored (in this case, FASTQ files).
    - The third volume specifies the directory where output files will be stored after the scripts have been executed. *You must create this directory before using the container to put files in it.*
* The `--user` argument identifies the user and group under which commands should be executed within the container. You should be able to leave this as is (if you are running on a Linux machine).
* The `--rm` argument indicates that Docker should automatically clean up the container and remove its file system when the container exits.
* `srp33/somatic_wgs` is the name (tag) of the Docker image; `latest` is a version tag associated with this image.
* The remaining arguments are specific to the task of using `bwa mem` to align FASTQ files to the reference genome.
    - The first argument (`-r`) indicates the name of a FASTA file that the user wishes to use as a reference genome.
    - The `-s1` and `-s2` arguments indicate the names of the FASTQ files that will be aligned; these files should be stored in the input_data volume specified above and should represent the first and second side of the paired-end reads, respectively.
    - The `-t` argument indicates the number of threads/cores that should be used during alignment; this argument should be a positive integer and is optional.
    - The `-o` argument indicates the name of the `bam` file that will be created (extension should be included).

## Feedback

[Let us know](https://github.com/srp33/docker_u54/issues) if you have any questions or problems.