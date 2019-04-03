#################### SVE IMAGE #######################
FROM wanpinglee/sve:0.1.0

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD TO PATH ########################
ENV PATH="/usr/local/bin/wgs:${PATH}"

################## INSTALL TOOLS ######################
ENV MINICONDA_VERSION=4.5.12
RUN apt-get update && apt-get install -y curl wget
RUN cd /
RUN curl -LO http://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh
RUN bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -p /miniconda -b
RUN rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh
ENV PATH=/miniconda/bin:${PATH}
RUN conda config --add channels bioconda \
 && conda config --add channels conda-forge
RUN conda install samtools=1.9 bwa=0.7.17 gatk4=4.1.0.0-0 delly=0.8.1-0 \
    sambamba=0.6.8 samblaster=0.1.24 manta=1.5.0 pysam=0.15.0 numpy=1.16.2
RUN conda create -n py2.7 python=2.7
RUN conda install strelka=2.9.10-0 lumpy-sv=0.2.13 subprocess32=3.5.3 htseq=0.11.2 cython=0.29.6 -n py2.7
ADD parliament2_overwrites/svtyper_env.yml /home/dnanexus/svtyper_env.yml
RUN conda env create -n svtyper_env --file /home/dnanexus/svtyper_env.yml

## TOOLS THAT DO NOT NEED TO BE INSTALLED UNDER PARLIAMENT2 #####
# py3.6 samtools sambamba samblaster
# py2.7 manta (maybe)

################## ADD SCRIPTS ########################
ADD echo_usage.sh /usr/local/bin/echo_usage
ADD usage_functions.sh /usr/local/bin/usage_functions
ADD check_functions.sh /usr/local/bin/check_functions
ADD call_somatic_variants_varscan.sh /usr/local/bin/call_somatic_variants_varscan
ADD check_permissions.py /

################# ADD PARLIAMENT2 OVERWRITES ##########
#ADD parliament2_overwrites/parliament2.py /home/dnanexus/parliament2.py
#ADD parliament2_overwrites/parliament2.sh /home/dnanexus/parliament2.sh
#ADD parliament2_overwrites/get_reference.py /home/dnanexus/get_reference.py
#ADD parliament2_overwrites/parallelize_svtyper.sh /home/dnanexus/parallelize_svtyper.sh

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
ADD run_fusorsv.sh /usr/local/bin/run_fusorsv
ADD run_sve.sh /usr/local/bin/run_sve
#RUN chmod a+rwx -R /home/dnanexus

################## SETUP WORKDIR #######################
WORKDIR /data

######## OVERRIDE PARLIAMENT2 ENTRYPOINT ##############
ENTRYPOINT [""]

##################### RUN BWA #########################
CMD echo_usage
