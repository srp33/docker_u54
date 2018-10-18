# Run bwa
#################### BWA IMAGE ########################
FROM biocontainers/biocontainers:latest

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## SET VARIABLES ######################
ENV THREADS=1
ENV REF_GENOME=Null
ENV SAMPLE=Null

################## ADD SCRIPT #########################
ADD run.sh /data/
ADD check_permissions.py /data/

################## INSTALL TOOLS ######################
RUN conda install bwa
RUN conda install samtools

################## SETUP WORKDIR #######################
WORKDIR /data

##################### RUN BWA #########################
CMD ./run.sh $REF_GENOME $THREADS $SAMPLE
