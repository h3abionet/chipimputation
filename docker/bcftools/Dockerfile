############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintenair="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ARG DEBIAN_FRONTEND=noninteractive

################## BEGIN INSTALLATION ######################

# Install wget
RUN apt-get update && apt-get install -y \
  autoconf \
  build-essential \
  git \
  libncurses5-dev \
  pkg-config \
  unzip \
  wget curl \
  python python-dev \
  libbz2-dev \
  liblzma-dev \
  zlib1g-dev &&\
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install htslib
RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
  tar -xvf htslib-1.9.tar.bz2 && \
  cd htslib-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf htslib-1.9*

# Install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
  tar -xvf samtools-1.9.tar.bz2 && \
  cd samtools-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf samtools-1.9*

# Install VCFTools
RUN wget https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz && \
  tar -xvf vcftools-0.1.16.tar.gz && \
  cd vcftools-0.1.16 && \
  ./configure && \
  make && \
  make install

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
  tar -xvf bcftools-1.13.tar.bz2 && \
  cd bcftools-1.13 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf bcftools-1.13*

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
  tar -zxvf bedtools-2.28.0.tar.gz && \
  cd bedtools2 && \
  make && \
  mv bin/bedtools /usr/local/bin/ && \
  cd .. && rm -r bedtools2

RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]