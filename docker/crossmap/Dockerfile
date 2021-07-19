############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 20.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintenair="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8

################## BEGIN INSTALLATION ######################

# Install wget
RUN apt-get update &&  \
  apt-get install -y python3.8 python3-pip && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN pip3 install CrossMap

RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
  tar -xvf htslib-1.9.tar.bz2 && \
  cd htslib-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf htslib-1.9*

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
  tar -xvf bcftools-1.13.tar.bz2 && \
  cd bcftools-1.13 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf bcftools-1.13*

USER ubuntu

CMD ["/bin/bash","-i"]