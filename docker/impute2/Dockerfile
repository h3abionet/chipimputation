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
    wget curl \
    python python-dev \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    r-base && \
  apt-get clean && \
  rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Install bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -zxvf bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make

# Install IMPUTE2
RUN wget http://mathgen.stats.ox.ac.uk/impute/impute_v2.3.2_x86_64_static.tgz && \
  tar -zxvf impute_v2.3.2_x86_64_static.tgz && \
  mv impute_v2.3.2_x86_64_static/impute2 /usr/local/bin/impute2 && \
  mkdir /opt/impute2/example -p && \
  mv impute_v2.3.2_x86_64_static/Example/* /opt/impute2/example && \
  rm -rf impute_v2.3.2_x86_64_static impute_v2.3.2_x86_64_static.tgz

# Install VCFTools
### Build from GitHub
RUN git clone https://github.com/vcftools/vcftools.git  && \
  cd vcftools && \
  ./autogen.sh && \
  ./configure && \
  make -j && \
  make install && \
  cd .. && rm -rf vcftools

# install R packages
RUN R --slave -e 'install.packages("dplyr", repos="https://cloud.r-project.org/")' && \
  R --slave -e 'install.packages("ggplot2", repos="https://cloud.r-project.org/")' && \
  R --slave -e 'install.packages("data.table", repos="https://cloud.r-project.org/")' && \
  R --slave -e 'install.packages("sm", repos="https://cloud.r-project.org/")' && \
  R --slave -e 'install.packages("optparse", repos="https://cloud.r-project.org/")'  && \
  R --slave -e 'install.packages("ggsci", repos="https://cloud.r-project.org/")'  && \
  R --slave -e 'install.packages("tidyr", repos="https://cloud.r-project.org/")'

# Install minimac3
RUN wget ftp://share.sph.umich.edu/minimac3/Minimac3.v2.0.1.tar.gz  && \
  tar -xzvf Minimac3.v2.0.1.tar.gz  && \
  cd Minimac3/  && \
  make  && \
  mv ./bin/Minimac3 /usr/local/bin/

# Install minimac4
RUN wget http://debian.mirror.ac.za/debian/pool/main/libs/libstatgen/libstatgen0_1.0.14-5_amd64.deb && \
  dpkg -i libstatgen0_1.0.14-5_amd64.deb && \
  wget http://debian.mirror.ac.za/debian/pool/main/m/minimac4/minimac4_1.0.0-2_amd64.deb && \
  dpkg -i minimac4_1.0.0-2_amd64.deb

# Install PLINK2
# there is an undocumented stable url (without the date)
RUN wget http://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_latest.zip -O plink.zip && \
  unzip plink.zip -d /usr/local/bin/ && \
  rm -f plink.zip

# Install bcftools
RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
  tar -xvf bcftools-1.9.tar.bz2 && \
  cd bcftools-1.9 && \
  ./configure --prefix=/usr/local && \
  make && \
  make install && \
  cd .. && rm -rf bcftools-1.9*

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

# Install Eagle
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz && \
    gunzip Eagle_v2.4.1.tar.gz && \
    tar xvf Eagle_v2.4.1.tar && \
    mv Eagle_v2.4.1/eagle /usr/local/bin/ && \
    rm -rf Eagle_v2.4.1


RUN useradd --create-home --shell /bin/bash ubuntu && \
  chown -R ubuntu:ubuntu /home/ubuntu

USER ubuntu

CMD ["/bin/bash","-i"]
