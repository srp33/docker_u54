# Run bwa
#################### BWA IMAGE ########################
FROM biocontainers/biocontainers:latest

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## SET VARIABLES ######################
ENV BWA_THREADS=100
ENV REF_GENOME=ucsc.hg19.fasta.gz

################## ADD SCRIPT #########################
ADD run.sh /data/

################## INSTALL TOOLS ######################
RUN conda install bwa
RUN conda install samtools

################## SETUP WORKDIR #######################
WORKDIR /data/results
WORKDIR /data

##################### RUN BWA #########################
CMD ./run.sh $REF_GENOME $BWA_THREADS
