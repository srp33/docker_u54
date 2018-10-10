# Run bwa
#################### BWA IMAGE ########################
FROM biocontainers/biocontainers:latest

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## SET VARIABLES ######################
ENV BWA_THREADS=100
ENV REF_GENOME=ucsc.hg19.fasta.gz

################## INSTALL TOOLS ######################
RUN conda install bwa
RUN conda install samtools

##################### RUN BWA #########################
CMD bwa mem /data/ref_index/$REF_GENOME /data/SampleData/* | samtools view -S -b
