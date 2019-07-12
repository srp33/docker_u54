FROM ubuntu:18.04

########### MODIFY ENVIRONMENT VARIABLES #############
ENV PATH="/usr/local/bin/wgs:/miniconda/bin:${PATH}"
ENV TZ=US
ENV MINICONDA_VERSION=4.6.14

################## INSTALL TOOLS ######################
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone \
  && apt-get update \
  && apt-get install -y curl wget bzip2 git-all build-essential zlib1g-dev locales \
  && locale-gen en_US.UTF-8 \
  && cd / \
  && curl -LO http://repo.continuum.io/miniconda/Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
  && bash Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh -p /miniconda -b \
  && rm Miniconda3-${MINICONDA_VERSION}-Linux-x86_64.sh \
  && conda config --add channels bioconda \
  && conda config --add channels conda-forge \
  && conda install samtools=1.9 bwa=0.7.17 gatk=3.8 gatk4=4.1.2.0-0 delly=0.8.1-0 \
                   sambamba=0.7.0 samblaster=0.1.24 pysam=0.15.0 numpy=1.16.2 \
                   picard=2.18.29-0 perl-list-moreutils=0.428 atropos=1.1.22-0 \
  && wget -O "GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2" "https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.8-1-0-gf15c1c3ef" \
  && tar -jxvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
  && rm GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 \
  && gatk3-register $(pwd)/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  && conda create -n py2.7 python=2.7 \
  && conda install strelka=2.9.10-0 lumpy-sv=0.2.13 subprocess32=3.5.3 \
                   htseq=0.11.2 cython=0.29.6 svtyper=0.7.0 survivor=1.0.6 manta=1.5.0 -n py2.7 \
  && alias awk=gawk \
  && mkdir -p /data \
  && chmod 777 /data

################## ADD SCRIPTS ########################
ADD src/*.py /
ADD src/system_commands/* /usr/local/bin/
ADD src/wgs/* /usr/local/bin/wgs/

################## SETUP WORKDIR #######################
WORKDIR /data

CMD echo_usage
