############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:20.04

# File Author / Maintainer
LABEL maintenair="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH
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

## Install IMPUTE5
RUN wget https://www.dropbox.com/sh/mwnceyhir8yze2j/AAD2VrkZze6ZLrcGX-jok4KRa/impute5_v1.1.5.zip && \
  unzip impute5_v1.1.5.zip && \
  mv impute5_v1.1.5/imp5Chunker_1.1.5_static /usr/local/bin/impute5Chunker && \
  mv impute5_v1.1.5/imp5Converter_1.1.5_static /usr/local/bin/impute5Converter && \
  mv impute5_v1.1.5/impute5_1.1.5_static /usr/local/bin/impute5 && \
  chmod a+x /usr/local/bin/imp5Chunker_1.1.5_static /usr/local/bin/imp5Converter_1.1.5_static /usr/local/bin/impute5_1.1.5_static && \
  rm -rf impute5_v1.1.5*

RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]
