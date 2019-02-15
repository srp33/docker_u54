# Run bwa
#################### GATK IMAGE #######################
FROM dnanexus/parliament2:v0.1.9-12-ge134fcda

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD TO PATH ########################
ENV PATH="/usr/local/bin/wgs:${PATH}"

################## INSTALL TOOLS ######################
#RUN conda config --add channels bioconda
RUN conda install bwa=0.7.17 gatk4=4.1.0.0-0 delly=0.8.1-0
RUN conda create -n py2.7 python=2.7
RUN conda install strelka=2.9.10-0 lumpy-sv=0.2.13 -n py2.7

## TOOLS THAT DO NOT NEED TO BE INSTALLED UNDER PARLIAMENT2 #####
# py3.6 samtools sambamba samblaster
# py2.7 manta (maybe)

################## ADD SCRIPTS ########################
ADD echo_usage.sh /usr/local/bin/echo_usage
ADD usage_functions.sh /usr/local/bin/usage_functions
ADD check_functions.sh /usr/local/bin/check_functions
ADD call_somatic_variants_varscan.sh /usr/local/bin/call_somatic_variants_varscan
ADD parliament2.py /home/dnanexus/parliament2.py
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
ADD base_recalibrator.sh /usr/local/bin/wgs/base_recalibrator
ADD apply_bqsr.sh /usr/local/bin/wgs/apply_bqsr
ADD call_structural_variants_manta.sh /usr/local/bin/wgs/call_structural_variants_manta
ADD samblast.sh /usr/local/bin/wgs/samblast
ADD call_structural_variants_lumpy.sh /usr/local/bin/wgs/call_structural_variants_lumpy
ADD call_structural_variants_delly.sh /usr/local/bin/wgs/call_structural_variants_delly

################ ADD OTHER SCRIPTS ####################
ADD run_parliament2.sh /usr/local/bin/run_parliament2

################## SETUP WORKDIR #######################
WORKDIR /data

######## OVERRIDE PARLIAMENT2 ENTRYPOINT ##############
ENTRYPOINT [""]

##################### RUN BWA #########################
CMD echo_usage
