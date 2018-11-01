# Run bwa
#################### BWA IMAGE ########################
FROM broadinstitute/gatk:latest

#################### MAINTAINER #######################
MAINTAINER Zachary Elias Ence <zac.ence@gmail.com>

################## ADD SCRIPT #########################
ADD align.sh /usr/local/bin/align
ADD echo_usage.sh /usr/local/bin/echo_usage
ADD check_permissions.py /
ADD echo_tester.sh /usr/local/bin/echo_tester

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
CMD ./run.sh $REF_GENOME $THREADS $SAMPLE
