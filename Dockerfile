# docker build --tag=bowtie-aligner .

FROM ubuntu:18.04
RUN apt-get update && apt-get install -y \
    curl \
    unzip \
    perl \
    default-jdk \
    wget \
    fastqc \
    samtools \
    sambamba \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python \
    python3.6 \
    python3-pip

######### Cutadapt Setup #############
RUN pip3 install --upgrade pip setuptools
RUN pip3 install --user --upgrade 'cutadapt==2.7'
RUN ln -s ~/.local/bin/cutadapt /usr/bin/
######### End Cutadapt Setup #########

######## Bowtie Setup ########
ARG APP_DIR=/usr/local/bin/
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie/1.2.2/bowtie-1.2.2-linux-x86_64.zip/download && unzip download && mv bowtie-1.2.2-linux-x86_64 ${APP_DIR}bowtie && rm download
ENV PATH ${APP_DIR}bowtie:$PATH
######## End Bowtie Setup ########

######### Pysam Setup ################
RUN pip3 install pysam
######### End Pysam Setup ############

######### UMI Tools Setup ############
RUN pip3 install --user --upgrade umi_tools
RUN ln -s ~/.local/bin/umi_tools /usr/bin/umi_tools
######### END UMI Tools Setup ########

######### Subread Setup #############
ENV SUBREAD_VERSION 1.6.4 
RUN mkdir -p /opt/subread
RUN curl -SLO https://sourceforge.net/projects/subread/files/subread-${SUBREAD_VERSION}/subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz \
    && tar -zxvf subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz --directory /opt/ \
    && rm subread-${SUBREAD_VERSION}-Linux-x86_64.tar.gz
ENV PATH /opt/subread-${SUBREAD_VERSION}-Linux-x86_64/bin/:$PATH
######### End Subread Setup ##########

######### R Setup ###############
ENV DEBIAN_FRONTEND=noninteractive 

RUN apt-get update && apt-get install -y \
    pandoc \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    gdebi-core -y

ENV R_VERSION=4.0.2
RUN curl -O https://cdn.rstudio.com/r/ubuntu-1804/pkgs/r-${R_VERSION}_1_amd64.deb
RUN gdebi -n r-${R_VERSION}_1_amd64.deb
RUN ln -s /opt/R/${R_VERSION}/bin/R /usr/local/bin/R
RUN ln -s /opt/R/${R_VERSION}/bin/Rscript /usr/local/bin/Rscript
RUN chmod +x /opt/R/${R_VERSION}/bin/R
RUN chmod +x /opt/R/${R_VERSION}/bin/Rscript
RUN Rscript -e 'install.packages(c("readr","dplyr","tidyr","ggplot2", "fastqcr","stringr","kableExtra","rlist","data.table","tibble"), repos = "http://cran.r-project.org")'
######### End R Setup ###########

######### Additional Python Setup ############
RUN pip3 install pandas --user
RUN pip3 install xlsxwriter --user
######### END Additional Python Setup ########

WORKDIR /opt/biorad

COPY . .
