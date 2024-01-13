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
    libfontconfig1-dev \
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

ENV urls="'https://cran.r-project.org/src/contrib/Archive/rlang/rlang_0.4.8.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/ellipsis/ellipsis_0.3.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/digest/digest_0.6.27.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/glue/glue_1.4.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/vctrs/vctrs_0.3.5.tar.gz',\
'https://cran.r-project.org/src/contrib/pkgconfig_2.0.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/hms/hms_0.5.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/lifecycle/lifecycle_0.2.0.tar.gz',\
'https://cran.r-project.org/src/contrib/assertthat_0.2.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/crayon/crayon_1.3.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/fansi/fansi_0.4.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/cli/cli_2.2.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/utf8/utf8_1.1.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/pillar/pillar_1.4.7.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/magrittr/magrittr_2.0.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/tibble/tibble_3.0.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/clipr/clipr_0.7.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/R6/R6_2.5.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/BH/BH_1.72.0-3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/cpp11/cpp11_0.2.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/readr/readr_1.4.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/purrr/purrr_0.3.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/tidyselect/tidyselect_1.1.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/ps/ps_1.4.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/processx/processx_3.4.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/callr/callr_3.5.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/rprojroot/rprojroot_2.0.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/desc/desc_1.2.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/withr/withr_2.3.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/prettyunits/prettyunits_1.1.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/pkgbuild/pkgbuild_1.1.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/rstudioapi/rstudioapi_0.13.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/pkgload/pkgload_1.1.0.tar.gz',\
'https://cran.r-project.org/src/contrib/rematch2_2.1.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/diffobj/diffobj_0.3.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/waldo/waldo_0.2.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/brio/brio_1.1.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/evaluate/evaluate_0.14.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/jsonlite/jsonlite_1.7.1.tar.gz',\
'https://cran.r-project.org/src/contrib/praise_1.0.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/testthat/testthat_3.0.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/isoband/isoband_0.2.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/gtable/gtable_0.3.0.tar.gz',\
'https://cran.r-project.org/src/contrib/gridExtra_2.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/farver/farver_2.0.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/labeling/labeling_0.4.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/colorspace/colorspace_2.0-0.tar.gz',\
'https://cran.r-project.org/src/contrib/munsell_0.5.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/RColorBrewer/RColorBrewer_1.1-2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/viridisLite/viridisLite_0.3.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/scales/scales_1.1.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/ggplot2/ggplot2_3.3.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/xfun/xfun_0.39.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/mime/mime_0.9.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/markdown/markdown_1.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/highr/highr_0.8.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/stringi/stringi_1.5.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/stringr/stringr_1.4.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/yaml/yaml_2.2.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/knitr/knitr_1.30.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/sys/sys_3.4.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/askpass/askpass_1.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/openssl/openssl_1.4.3.tar.gz',\
'https://cran.r-project.org/src/contrib/selectr_0.4-2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/xml2/xml2_1.3.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/curl/curl_4.3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/httr/httr_1.4.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/rvest/rvest_0.3.6.tar.gz',\
'https://cran.r-project.org/src/contrib/base64enc_0.1-3.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/htmltools/htmltools_0.5.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/tinytex/tinytex_0.27.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/rmarkdown/rmarkdown_2.5.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/generics/generics_0.1.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/dplyr/dplyr_1.0.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/tidyr/tidyr_1.1.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/webshot/webshot_0.5.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/KernSmooth/KernSmooth_2.23-17.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-51.6.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.2-18.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/XML/XML_3.99-0.5.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/boot/boot_1.3-25.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/class/class_7.3-17.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/cluster/cluster_2.1.0.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/codetools/codetools_0.2-16.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/data.table/data.table_1.13.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/fastqcr/fastqcr_0.1.2.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/foreign/foreign_0.8-80.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/kableExtra/kableExtra_1.3.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/lattice/lattice_0.20-41.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/mgcv/mgcv_1.8-31.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/nlme/nlme_3.1-148.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/nnet/nnet_7.3-14.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/rlist/rlist_0.4.6.1.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/rpart/rpart_4.1-15.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/spatial/spatial_7.3-12.tar.gz',\
'https://cran.r-project.org/src/contrib/Archive/survival/survival_3.1-12.tar.gz'"
RUN Rscript -e "install.packages(c($urls),repos=NULL,type='source')"

#RUN Rscript -e 'install.packages(c(),repos=NULL, type="source")'
#RUN Rscript -e 'install.packages(c("readr","dplyr","tidyr","ggplot2", "fastqcr","stringr","kableExtra","rlist","data.table","tibble"), repos = "http://cran.r-project.org")'
######### End R Setup ###########

######### Additional Python Setup ############
RUN pip3 install pandas --user
RUN pip3 install xlsxwriter --user
######### END Additional Python Setup ########

WORKDIR /opt/biorad

COPY . .
