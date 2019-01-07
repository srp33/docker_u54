# Run bwa
#################### GATK IMAGE ########################
FROM broadinstitute/gatk:4.0.11.0

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD SCRIPTS ########################
ADD echo_usage.sh /usr/local/bin/echo_usage
ADD usage_functions.sh /usr/local/bin/usage_functions
ADD check_functions.sh /usr/local/bin/check_functions
ADD call_somatic_variants_varscan.sh /usr/local/bin/call_somatic_variants_varscan
ADD check_permissions.py /

################ ADD WGS SCRIPTS ######################
ADD bwa_mem_align.sh /usr/local/bin/wgs/bwa_mem_align
ADD call_somatic_variants_gatk4.sh /usr/local/bin/wgs/call_somatic_variants_gatk4
ADD call_somatic_variants_strelka.sh /usr/local/bin/wgs/call_somatic_variants_strelka
ADD index_bam.sh /usr/local/bin/wgs/index_bam
ADD mark_duplicates.sh /usr/local/bin/wgs/mark_duplicates
ADD merge_bams.sh /usr/local/bin/wgs/merge_bams
ADD slice_bam.sh /usr/local/bin/wgs/slice_bam
ADD sort_bam.sh /usr/local/bin/wgs/sort_bam
ADD add_read_groups.sh /usr/local/bin/wgs/add_read_groups
ADD samtools_mpileup.sh /usr/local/bin/wgs/samtools_mpileup

################## ADD TO PATH ########################
ENV PATH="/usr/local/bin/wgs:${PATH}"

################## INSTALL TOOLS ######################
RUN conda config --add channels bioconda
RUN conda install bwa samtools sambamba varscan
RUN conda install strelka
#RUN mkdir -p /tmp && chmod 755 /data

################## SETUP WORKDIR #######################
WORKDIR /data

##################### RUN BWA #########################
CMD echo_usage
