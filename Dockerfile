FROM ubuntu:22.04


ENV DEBIAN_FRONTEND=noninteractive


RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gnupg \
    apt-transport-https \
    ca-certificates \
    python3-pip \
    python3-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libssl-dev \
    && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
    && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' \
    && apt-get update && apt-get install -y r-base r-base-dev


RUN pip3 install --upgrade pip setuptools wheel
RUN pip3 install pandas pytest rpy2==3.5.1


RUN R -e "install.packages('BiocManager', repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install(c('DESeq2', 'airway'))"


ENV R_HOME=/usr/lib/R
ENV LD_LIBRARY_PATH=/usr/lib/R/lib:/usr/local/lib:/usr/lib/x86_64-linux-gnu

WORKDIR /app
CMD ["/bin/bash"]