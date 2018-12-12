# Run bwa
#################### GATK IMAGE ########################
FROM broadinstitute/gatk:4.0.11.0

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD SCRIPT #########################
ADD bwa_mem_align.sh /usr/local/bin/bwa_mem_align
ADD echo_usage.sh /usr/local/bin/echo_usage
ADD usage_functions.sh /usr/local/bin/usage_functions
ADD check_functions.sh /usr/local/bin/check_functions
ADD call_gatk_variants.sh /usr/local/bin/call_gatk_variants
ADD call_somatic_variants_strelka.sh /usr/local/bin/call_somatic_variants_strelka
ADD call_somatic_variants_varscan.sh /usr/local/bin/call_somatic_variants_varscan
ADD index_bam.sh /usr/local/bin/index_bam
ADD mark_duplicates.sh /usr/local/bin/mark_duplicates
ADD merge_bams.sh /usr/local/bin/merge_bams
ADD slice_bam.sh /usr/local/bin/slice_bam
ADD sort_bam.sh /usr/local/bin/sort_bam
ADD add_read_groups.sh /usr/local/bin/add_read_groups
ADD samtools_mpileup.sh /usr/local/bin/samtools_mpileup
ADD check_permissions.py /

################## INSTALL TOOLS ######################
RUN conda config --add channels bioconda
RUN conda install bwa samtools sambamba varscan
RUN conda install strelka
#RUN mkdir -p /tmp && chmod 755 /data

################## SETUP WORKDIR #######################
WORKDIR /data

##################### RUN BWA #########################
CMD echo_usage
