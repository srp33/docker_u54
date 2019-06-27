FROM ubuntu:18.04

########### MODIFY ENVIRONMENT VARIABLES #############
ENV PATH="/usr/local/bin/wgs:/miniconda/bin:${PATH}"
ENV TZ=US
ENV MINICONDA_VERSION=4.5.12

################## INSTALL TOOLS ######################
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone \
  && apt-get update \
  && apt-get install -y curl wget bzip2 git-all build-essential zlib1g-dev \
  && cd / \
  && curl -LO http://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
  && bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -p /miniconda -b \
  && rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge \
  && conda install samtools=1.9 bwa=0.7.17 gatk4=4.1.0.0-0 delly=0.8.1-0 \
                   sambamba=0.6.8 samblaster=0.1.24 manta=1.5.0 pysam=0.15.0 numpy=1.16.2 \
                   picard=2.18.29-0 perl-list-moreutils=0.428 \
  && conda create -n py2.7 python=2.7 \
  && conda install strelka=2.9.10-0 lumpy-sv=0.2.13 subprocess32=3.5.3 \
                   htseq=0.11.2 cython=0.29.6 svtyper=0.7.0 survivor=1.0.6 -n py2.7 \
  && alias awk=gawk

################## ADD SCRIPTS ########################
ADD src/system_commands/echo_usage.sh /usr/local/bin/echo_usage
ADD src/system_commands/usage_functions.sh /usr/local/bin/usage_functions
ADD src/system_commands/check_functions.sh /usr/local/bin/check_functions
ADD src/system_commands/entrypoint.sh /usr/local/bin/wgs/entrypoint
ADD src/system_commands/check_permissions.py /

################ ADD WGS SCRIPTS ######################
ADD src/wgs/bwa_mem_align.sh /usr/local/bin/wgs/bwa_mem_align
ADD src/wgs/call_somatic_variants_gatk4.sh /usr/local/bin/wgs/call_somatic_variants_gatk4
ADD src/wgs/call_somatic_variants_strelka.sh /usr/local/bin/wgs/call_somatic_variants_strelka
ADD src/wgs/index_bam.sh /usr/local/bin/wgs/index_bam
ADD src/wgs/mark_duplicates.sh /usr/local/bin/wgs/mark_duplicates
ADD src/wgs/merge_bams.sh /usr/local/bin/wgs/merge_bams
ADD src/wgs/slice_bam.sh /usr/local/bin/wgs/slice_bam
ADD src/wgs/sort_bam.sh /usr/local/bin/wgs/sort_bam
ADD src/wgs/add_read_groups.sh /usr/local/bin/wgs/add_read_groups
ADD src/wgs/samtools_mpileup.sh /usr/local/bin/wgs/samtools_mpileup
ADD src/wgs/base_recalibrator.sh /usr/local/bin/wgs/base_recalibrator
ADD src/wgs/apply_bqsr.sh /usr/local/bin/wgs/apply_bqsr
ADD src/wgs/call_structural_variants_manta.sh /usr/local/bin/wgs/call_structural_variants_manta
ADD src/wgs/samblast.sh /usr/local/bin/wgs/samblast
ADD src/wgs/call_structural_variants_lumpy.sh /usr/local/bin/wgs/call_structural_variants_lumpy
ADD src/wgs/call_structural_variants_delly.sh /usr/local/bin/wgs/call_structural_variants_delly
ADD src/wgs/svtype_vcf.sh /usr/local/bin/wgs/svtype_vcf
ADD src/wgs/run_survivor.sh /usr/local/bin/wgs/run_survivor

################## SETUP WORKDIR #######################
WORKDIR /data
RUN chmod 777 /data

CMD echo_usage
