#! /bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This command is not available to users" && exit 1
fi

usage_align (){
echo "bwa_mem_align

Description:
Align FASTQ files to a reference genome using the Burrows-Wheeler Aligner software.

Options:
  -r, --reference <name of reference genome FASTA file>
  -s1, --sample1 <file 1>
  -s2, --sample2 <file 2>
  -o, --output <name of outputted BAM file>
  -h, --help
  -t, --nthreads <number of threads> (Optional)

Usage:
docker run \\
  -v <location of reference FASTA file>:/data/ref_genome \\
  -v <location of FASTQ files>:/data/input_data \\
  -v <location for outputted BAM file>:/data/output_data \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  bwa_mem_align \\
    -r <reference FASTA file> \\
    -s1 <file 1> \\
    -s2 <file 2> \\
    -o <name of outputted BAM file> \\
    -t <number of threads> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of reference FASTA file>
  <location of FASTQ files>
  <location for outputted BAM file>
"
}

usage_index_bam (){
echo "index_bam

Description:
Index a BAM file.

Options:
  -t, --nthreads <number of threads> (Optional)
  -b, --bam <name of BAM file>
  -h, --help

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  index_bam \\
     -b <BAM file> \\
     -t <number of threads> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
"
}

usage_mark_duplicates (){
echo "mark_duplicates

Description:
Mark duplicates in a BAM file.

Options:
  -t, --nthreads <number of threads> (Optional)
  -b, --bam <name of BAM file>
  -o, --output <name of output file>
  -h, --help

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  mark_duplicates \\
    -b <BAM file> \\
    -o <name of output file> \\
    -t <number of threads>

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}

usage_merge_bams (){
echo "merge_bams

Description:
Merge multiple BAM files into one and index it

Options:
  -t, --nthreads <number of threads> (Optional)
  -b, --bam <name of BAM file> (2 or more required)
  -o, --output <name of output file>
  -h, --help

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  merge_bams \\
    -b <BAM file> \\
    -b <BAM file> \\
    ... \\
    -o <name of output file> \\
    -t <number of threads>

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}

usage_slice_bam (){
echo "slice_bam

Description:
Slice a BAM file based on chromosomal coordinates.

Options:
  -t, --nthreads <number of threads> (Optional)
  -b, --bam <name of BAM file>
  -r, --region <region to slice> (e.g. \"chr2\" or \"chr2:1000-2000\")
  -o, --output <name of output file>
  -h, --help

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  slice_bam \\
    -b <BAM file> \\
    -r <region to slice> \\
    -o <name of output file> \\
    -t <number of threads>

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}

usage_sort_bam (){
echo "sort_bam

Description:
Sort and index a BAM file.

Options:
  -t, --nthreads <number of threads> (Optional)
  -b, --bam <name of BAM file>
  -o, --output <name of output file>
  -h, --help

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  sort_bam \\
    -b <BAM file> \\
    -o <name of output file> \\
    -t <number of threads>

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}

usage_strelka (){
echo "call_somatic_variants_strelka

Options:
  -b, --bam <name of BAM file>
  -r, --reference <name of reference genome FASTA file>
  -h, --help

ERROR: THIS COMMAND IS NOT READY FOR USAGE
"
}

usage_varscan (){
echo "call_somatic_variants_varscan

Options:
  -p --pileup <name of pileup file>
  -o, --output <name of output file>
  -h, --help

ERROR: THIS COMMAND IS NOT READY FOR USAGE
"
}

usage_mutect (){
echo "call_gatk_variants

Options:
  -t, --tumor <name of tumor BAM file>
  -n, --normal <name of normal BAM file
  -o, --output <name of output file>
  -r, --reference <name of reference file>
  -h, --help

ERROR: THIS COMMAND IS NOT READY FOR USAGE
"
}

usage_add_read_groups (){
echo "add_read_groups

Description:
Add/replace read group in BAM file header

Options:
  -b, --bam_file <name of BAM file>
  -id, --group_id <group_id>
  -lb, --library_identifier <name of library>
  -s, --sample <name of sample>
  -o, --output <name of output file>
  -h, --help

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  ---user $(id -u):$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  add_read_groups \\
    -b <BAM file> \\
    -id <Group ID> \\
    -lb <Library Identifier> \\
    -s <Sample> \\
    -o <name of output file>

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}