############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 20.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:20.04

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

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
  /bin/bash ~/miniconda.sh -b -p /opt/conda && \
  rm ~/miniconda.sh && \
  echo "conda activate base" >> ~/.bashrc

RUN conda update conda && \
  conda clean --all --yes && \
  conda install -y -c bioconda crossmap


RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]