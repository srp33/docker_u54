#! /bin/bash

if [[ "${BASH_SOURCE[0]}" == "${0}" ]]; then
    echo "ERROR: This command is not available to users" && exit 1
fi

usage_align (){
echo "
bwa_mem_align
Options:
-t, --nthreads <number of threads> (Optional)
-r, --reference <name of reference genome fasta file>
-s1, --sample1 <read 1>
-s2, --sample2 <read 2>
-h, --help
"
}

usage_index_bam (){
echo "
index_bam
Options:
-t, --nthreads <number of threads> (Optional)
-b, --bam <name of bam file>
-h, --help
"
}

usage_mark_duplicates (){
echo "
mark_duplicates
Options:
-t, --nthreads <number of threads> (Optional)
-b, --bam <name of bam file>
-o, --output <name of output file>
-h, --help
"
}

usage_merge_bams (){
echo "
merge_bams
Options:
-t, --nthreads <number of threads> (Optional)
-b, --bam <name of bam file>
-o, --output <name of output file>
-h, --help
"
}

usage_slice_bam (){
echo "
slice_bam
Options:
-t, --nthreads <number of threads> (Optional)
-b, --bam <name of bam file>
-r, --region <region to slice>
-o, --output <name of output file>
-h, --help
"
}

usage_sort_bam (){
echo "
sort_bam
Options:
-t, --nthreads <number of threads> (Optional)
-b, --bam <name of bam file>
-h, --help
"
}