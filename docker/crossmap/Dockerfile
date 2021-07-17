############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
LABEL maintenair="Mamana Mbiyavanga mamana.mbiyavanga@uct.ac.za"

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

################## BEGIN INSTALLATION ######################

# Install wget
RUN apt-get update && apt-get install -y \
  autoconf \
  build-essential \
  git \
  pkg-config \
  unzip \
  wget curl \
  python python-dev &&\
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN conda update conda && \
  conda clean --all --yes && \
  conda install -y -c bioconda crossmap


RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]