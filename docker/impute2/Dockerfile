############################################################
# Dockerfile to build Genotype imputation
# Based on Ubuntu 16.04
############################################################

# Set the base image to Ubuntu
FROM ubuntu:16.04

# File Author / Maintainer
MAINTAINER Mamana Mbiyavanga "mamana.mbiyavanga@uct.ac.za"


################## BEGIN INSTALLATION ######################
# Install Basic tools

# Install wget
RUN apt-get update && apt-get install -y \
      autoconf \
      build-essential \
      git \
      libncurses5-dev \
      pkg-config \
      unzip \
      wget \
      python python-dev \
      zlib1g-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install IMPUTE2
RUN wget http://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
    tar -zxvf impute_v2.3.2_x86_64_static.tgz && \
        mv impute_v2.3.2_x86_64_static/impute2 /usr/local/bin/impute2 && \
           mkdir /opt/impute2/example -p && \
                 mv impute_v2.3.2_x86_64_static/Example/* /opt/impute2/example && \
                    rm -rf impute_v2.3.2_x86_64_static impute_v2.3.2_x86_64_static.tgz

# Install PLINK2
RUN wget https://www.cog-genomics.org/static/bin/plink180913/plink_linux_x86_64.zip -O plink.zip && \
    unzip plink.zip -d /usr/local/bin/ && \
    rm -f plink.zip

# Install VCFTools
### Build from GitHub
RUN git clone https://github.com/vcftools/vcftools.git  && \
    cd vcftools && \
    ./autogen.sh && \
    ./configure && \
    make -j && \
    make install && \
    cd .. && rm -rf vcftools

# Install Eagle
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.tar.gz && \
    gunzip Eagle_v2.4.tar.gz && \
    tar xvf Eagle_v2.4.tar && \
    mv Eagle_v2.4/eagle /usr/local/bin/ && \
    rm -rf Eagle_v2.4

RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]
