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
  -s2, --sample2 <file 2> (Optional)
  -o, --output <name of outputted BAM file>
  -h, --help
  -t, --nthreads <number of threads> (Optional)
  -c, --nchunks <number of chunks> (Optional)
  -p, --process_chunk <chunk to be processed> (Optional)
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of reference FASTA file>:/data/ref_genome \\
  -v <location of reference index files>:/data/ref_index \\
  -v <location of FASTQ files>:/data/input_data \\
  -v <location for outputted BAM file>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  bwa_mem_align \\
    -r <reference FASTA file> \\
    -s1 <file 1> \\
    -s2 <file 2> (Optional) \\
    -o <name of outputted BAM file> \\
    -t <number of threads> (Optional) \\
    -c <number of chunks> (Optional) \\
    -p <chunk to be processed (Optional) \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of reference FASTA file>
  <location of FASTQ files>
  <location for outputted BAM file>

  The --nchunks argument divides the FASTQ files into n number of chunks. The --process_chunk argument \
specifies which chunk should be aligned. --process_chunk must be less than --nchunks and greater than or \
equal to 0 (the default for --process_chunk is 0). --nchunks must be greater than 0, however, inputting \
1 for --nchunks will not change anything since that is the default behavior.
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
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for version log>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  index_bam \\
     -b <BAM file> \\
     -t <number of threads> (Optional)
     --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for version log>

  --log is an optional argument, meaning that the volume \
<location for version log>:/data/output_data is only required if --log is called.
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
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  mark_duplicates \\
    -b <BAM file> \\
    -o <name of output file> \\
    -t <number of threads> (Optional) \\
    --log <destination file for log> (Optional)

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
  -b, --bam <name of BAM file> (2 or more required)
  -o, --output <name of output file>
  -h, --help
  -t, --nthreads <number of threads> (Optional)
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  merge_bams \\
    -b <BAM file> \\
    -b <BAM file> \\
    ... \\
    -o <name of output file> \\
    -t <number of threads> (Optional) \\
    --log <destination file for log> (Optional)

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
  -b, --bam <name of BAM file>
  -r, --region <region to slice> (e.g. \"chr2\" or \"chr2:1000-2000\")
  -o, --output <name of output file>
  -h, --help
  -t, --nthreads <number of threads> (Optional)
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  slice_bam \\
    -b <BAM file> \\
    -r <region to slice> \\
    -o <name of output file> \\
    -t <number of threads> (Optional) \\
    --log <destination file for log> (Optional)

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
  -b, --bam <name of BAM file>
  -o, --output <name of output file>
  -h, --help
  -t, --nthreads <number of threads> (Optional)
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  sort_bam \\
    -b <BAM file> \\
    -o <name of output file> \\
    -t <number of threads> (Optional) \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}

usage_strelka (){
echo "call_somatic_variants_strelka

Description:
Calls somatic variants using Strelka for a tumor/normal pair

Options:
  -t, --tumorBam <name of tumor BAM file>
  -n, --normalBam <name of normal BAM file>
  -r, --referenceFasta <name of reference genome FASTA file>
  -i, --indelCandidates <name of VCF of indel alleles> (Optional)
  -c, --callRegions <name of file containing regions to call> (Optional)
  -d, --runDir <desired name for directory to be created where
                workflow scripts and output will be written> (Optional)
                [Default: StrelkaSomaticWorkflow]
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of indel alleles VCF and call regions file>:/data/input_data \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location of FASTA reference file>:/data/ref_genome \\
  -v <location of .fai reference index file>:/data/ref_index \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  call_somatic_variants_strelka \\
    -t <tumor BAM file> \\
    -n <normal BAM file> \\
    -r <reference FASTA file> \\
    -i <indel candidates VCF file> (Optional) \\
    -c <call regions file> (Optional) \\
    -d <output directory name> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of indel alleles VCF and call regions file>
  <location of BAM files>
  <location of FASTA reference file>
  <location of .fai reference index file>
  <location for output>

  This command currently requires tumor and normal BAM files.

  It should also be noted that strelka requires a .fai index file for the reference genome. \
If this file cannot be found in the /data/ref_index volume, it will be created. Unfortunately,\
 samtools cannot index gzipped files. Thus gzipped fasta files will be gunzipped into a \
temporary directory within the container.

  The indel candidates VCF file can be created using the call_structural_variants_manta \
command.
"
}

usage_manta (){
echo "call_structural_variants_manta

Description:
Calls somatic variants using Manta for a tumor/normal pair

Options:
  -t, --tumorBam <name of tumor BAM file>
  -n, --normalBam <name of normal BAM file>
  -r, --referenceFasta <name of reference genome FASTA file>
  -d, --runDir <desired name for directory to be created where
                workflow scripts and output will be written> (Optional)
                [Default: MantaWorkflow]
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location of FASTA reference file>:/data/ref_genome \\
  -v <location of .fai reference index file>:/data/ref_index \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  call_structural_variants_manta \\
    -t <tumor BAM file> \\
    -n <normal BAM file> \\
    -r <reference FASTA file> \\
    -d <output directory name> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location of FASTA reference file>
  <location of .fai reference index file>
  <location for output>

  This command currently requires tumor and normal BAM files. Both files must have been \
previously sorted. This can be done using the sort_bam command.

  It should also be noted that manta requires a .fai index file for the reference genome. \
If this file cannot be found in the /data/ref_index volume, it will be created. Unfortunately,\
 samtools cannot index gzipped files. Thus gzipped fasta files will be gunzipped into a \
temporary directory within the container.
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
echo "call_somatic_variants_gatk4

Description:
Calls somatic variants using Mutect2 for a tumor/normal pair

Options:
  -t, --tumor <name of tumor BAM file>
  -ts, --tumor_sample <name of tumor sample> (Optional)
  -n, --normal <name of normal BAM file>
  -ns --normal_sample <name of normal sample> (Optional)
  -o, --output <name of output file>
  -r, --reference <name of reference FASTA file>
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location of FASTA reference file>:/data/ref_genome \\
  -v <location of .fai reference index file>:/data/ref_index \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  call_somatic_variants_gatk4 \\
    -t <tumor BAM file> \\
    -ts <tumor sample> (Optional) \\
    -n <normal BAM file> \\
    -ns <normal sample> (Optional) \\
    -r <reference FASTA file> \\
    -o <name of output BAM file> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location of FASTA reference file>
  <location for output>

  Read group information must be included for in the tumor and normal BAM files. \
See 'add_read_groups' for details on how to add this information.

  It should also be noted that Mutect2 requires a .fai index file for the reference genome. \
If this file cannot be found in the /data/ref_index volume, it will be created. Unfortunately,\
 samtools cannot index gzipped files. Thus gzipped fasta files will be gunzipped into a \
temporary directory within the container.

  Tumor and normal samples should be provided, but if they are not, they will be taken from the names of \
the respective BAM files.
"
}

usage_add_read_groups (){
echo "add_read_groups

Description:
Add/replace read groups in BAM file

Options:
  -b, --bam_file <name of BAM file>
  -id, --group_id <group_id>
  -lb, --library_identifier <name of library>
  -s, --sample <name of sample>
  -o, --output <name of output file>
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  add_read_groups \\
    -b <BAM file> \\
    -id <Group ID> \\
    -lb <Library Identifier> \\
    -s <Sample> \\
    -o <name of output file> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location for output>
"
}

usage_mpileup () {
echo "add_read_groups

Description:
Create samtools pileup file

Options:
  -t, --tumor <name of tumor BAM file>
  -n, --normal <name of normal BAM file>
  -r, --reference <name of reference FASTA file>
  -o, --output <name of output file>
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location of reference FASTA file>:/data/ref_genome \\
  -v <location of reference .fai index file>:/data/ref_index \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  samtools_mpileup \\
    -t <tumor BAM file> \\
    -n <normal BAM file> \\
    -r <reference FASTA file> \\
    -o <name of output file> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location of reference FASTA file>
  <location of .fai index file>
  <location for output>

  If .fai is not found, it will be created. If reference FASTA file is gzipped, a temporary \
copy will be gunzipped into the container for the duration of this process. This will \
considerably lengthen the process.
"
}

usage_apply_bqsr () {
echo "apply_bqsr

Description:
Apply base quality score recalibration

Options:
  -b, --bam_file <name of input BAM file>
  -bqsr, --bqsr_recal_file <name of input recalibration table for BQSR>
  -o, --output <name of output BAM file>
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location of input recalibration table for BQSR>:/data/input_data \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  apply_bqsr \\
    -b <name of input BAM file> \\
    -bqsr <name of recalibration table> \\
    -o <name of output BAM file> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location of input recalibration table for BQSR>
  <location for output>
"
}

usage_base_recalibrator () {
echo "base_recalibrator

Description:
Detect systematic errors in base quality scores

Options:
  -b, --bam_file <name of input BAM file>
  -r, --reference <name of reference FASTA file>
  -s, --known_sites <URL to database of known polymorphic sites> (May be called more than once)
  -o, --output <name of output recalibration table file>
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM files>:/data/bam_files \\
  -v <location of VCF known sites files>:/data/vcf_files \\
  -v <location of reference FASTA file>:/data/ref_genome \\
  -v <location of reference .fai index file>:/data/ref_index \\
  -v <location for output>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  base_recalibrator \\
    -b <name of input BAM file \\
    -r <reference FASTA file> \\
    -s <name of VCF file> \\
    -s <name of VCF file> \\
    -o <name of output recalibration table file> \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM files>
  <location of VCF known sites files>
  <location of reference FASTA file>
  <location of .fai index file>
  <location for output>

  If .fai is not found, it will be created. If reference FASTA file is gzipped, a temporary \
copy will be gunzipped into the container for the duration of this process. This will \
considerably lengthen the process.

  If .tbi of VCF files are not found, they will be created.
"
}

usage_samblast () {
echo "samblast

Description:
Create split read and discordant read SAM files from BAM file

Options:
  -b, --bam_file <name of input BAM file>
  -o, --output <destination for output> (Default: /dev/null)
  -h, --help
  --log <destination file for log> (Optional)

Usage:
docker run \\
  -v <location of BAM file>:/data/bam_files \\
  -v <location for output and SAM files>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  samblast \\
    -b <name of input BAM file> \\
    -o <destination for output> (Default: /dev/null) \\
    --log <destination file for log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM file>
  <location for output and SAM files>
"
}

usage_lumpy () {
echo "call_structural_variants_lumpy

Description:
Call structural variants using lumpy

Options:
  -b, --bam_file <name of input BAM file>
  -s, --split_reads <name of split-reads BAM file>
  -d, --discordant_reads <name of discordant-reads BAM file>
  -o, --output <name of output VCF file>
  -h, --help
  --log <destination file for version log> (Optional)

Usage:
docker run \\
  -v <location of BAM file>:/data/bam_files \\
  -v <location for output VCF file>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  call_structural_variants_lumpy \\
    -b <name of input BAM file> \\
    -s <name of split-reads BAM file> \\
    -d <name -f discordant-reads BAM file> \\
    -o <name of output VCF file> \\
    --log <destination file for version log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM file>
  <location for output VCF file>

  All BAM files must be coordinate-sorted. If not, Lumpy will throw an error. Please ensure that all \
BAM files are sorted by using sort_bam before using call_structural_variants_lumpy.
"
}

usage_delly () {

    echo "call_structural_variants_delly

Description:
Call structural variants using delly

Options:
  -b, --bam_file <name of input BAM file>
  -r, --reference <name of reference FASTA file>
  -o, --output <name of output VCF file>
  -h, --help
  --log <destination file for version log> (Optional)

Usage:
docker run \\
  -v <location of BAM file>:/data/bam_files \\
  -v <location of reference FASTA file>:/data/ref_genome \\
  -v <location of reference .fai index file>:/data/ref_index \\
  -v <location for output VCF file>:/data/output_data \\
  --user \$(id -u):\$(id -g) \\
  --rm \\
  srp33/somatic_wgs:latest \\
  call_structural_variants_delly \\
    -b <name of input BAM file> \\
    -r <name of reference FASTA file> \\
    -o <name of output VCF file> \\
    --log <destination file for version log> (Optional)

Notes:

  To avoid permissions issues, please ensure that the following directories have been \
created on the host operating system before executing this command:

  <location of BAM file>
  <location of reference FASTA file>
  <location of reference .fai index file>
  <location for output VCF file>

  From the Delly docs: \"Delly needs a sorted, indexed and duplicate marked \
bam file for every input sample.\"
"
}