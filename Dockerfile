FROM rocker/r-ver:3.6.1

RUN apt-get update --fix-missing -qq && \
    apt-get install -y -q \
    vim \
    git \
    python3 \
    python3-pip \
    python \
    python-pip \
    libz-dev \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libpng-dev \
    libjpeg-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libgl-dev \
    libgsl-dev \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
# sys collections
RUN pip3 install numpy cython
RUN pip3 install HTSeq
#
RUN R -e 'install.packages(c("BiocManager", "devtools", "argparse", "dbplyr"))'
#
RUN R -e 'BiocManager::install("GenomicFeatures")'
#
ADD scripts /home/scripts
ADD software /home/software
RUN mkdir /home/data




