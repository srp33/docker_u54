# Run bwa
#################### BWA IMAGE ########################
FROM broadinstitute/gatk:latest

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD SCRIPT #########################
ADD bwa_mem_align.sh /usr/local/bin/bwa_mem_align
ADD echo_usage.sh /usr/local/bin/echo_usage
ADD usage_functions.sh /usr/local/bin/usage_functions
ADD check_for_args.sh /usr/local/bin/check_for_args
ADD call_gatk_variants.sh /usr/local/bin/call_gatk_variants
ADD call_somatic_variants_strelka.sh /usr/local/bin/call_somatic_variants_strelka
ADD call_somatic_variants_varscan.sh /usr/local/bin/call_somatic_variants_varscan
ADD index_bam.sh /usr/local/bin/index_bam
ADD mark_duplicates.sh /usr/local/bin/mark_duplicates
ADD merge_bams.sh /usr/local/bin/merge_bams
ADD slice_bam.sh /usr/local/bin/slice_bam
ADD sort_bam.sh /usr/local/bin/sort_bam
ADD strelka.sh /usr/local/bin/strelka
ADD check_permissions.py /

################## INSTALL TOOLS ######################
RUN conda config --add channels bioconda
RUN conda install bwa
RUN conda install samtools
RUN conda install sambamba
RUN conda install varscan
RUN conda install strelka

################## SETUP WORKDIR #######################
WORKDIR /data

##################### RUN BWA #########################
CMD echo_usage
