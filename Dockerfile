FROM rocker/r-ver:4.3.1

#  Install system dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential=12.9ubuntu3 \
    curl=7.81.0-1ubuntu1.20 \
    git \
    unzip=6.0-26ubuntu3.2 \
    libcurl4-openssl-dev=7.81.0-1ubuntu1.20 \
    libssl-dev=3.0.2-0ubuntu1.19  \
    libxml2-dev=2.9.13+dfsg-1ubuntu0.7 \
    libgit2-dev=1.1.0+dfsg.1-4.1ubuntu0.1 \
    libfreetype6-dev=2.11.1+dfsg-1ubuntu0.3 \
    libpng-dev=1.6.37-3build5 \
    libtiff5-dev=4.3.0-6ubuntu0.10 \
    libharfbuzz-dev=2.7.4-1ubuntu3.2 \
    libfribidi-dev=1.0.8-2ubuntu3.1 \
    libfontconfig1-dev=2.13.1-4.2ubuntu5 \
    cmake=3.22.1-1ubuntu1.22.04.2 \
    make=4.3-4.1build1 \
    g++=4:11.2.0-1ubuntu1 \
    zlib1g-dev=1:1.2.11.dfsg-2ubuntu9.2 \
    wget=1.21.2-2ubuntu1.1 \
    bedtools=2.30.0+dfsg-2ubuntu0.1 \
    samtools=1.13-4 \
    && apt-get clean

# Install LiftOver
RUN wget -q http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod +x liftOver && \
    mv liftOver /usr/local/bin/liftOver

# Install ichorCNA
RUN git clone https://github.com/broadinstitute/ichorCNA.git && \
    R -e "install.packages('plyr', version='1.8.9', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('BiocManager', version='1.30.26', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('remotes', version='2.5.0', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "install.packages('optparse', version='1.7.5', dependencies=TRUE, repos='http://cran.rstudio.com/')" && \
    R -e "BiocManager::install(c('HMMcopy','GenomeInfoDB','GenomicRanges'), ask=FALSE)" && \
    R -e "remotes::install_local('ichorCNA', dependencies=TRUE)"

# Download UCSC liftover chain file and chrom sizes for hg38
RUN mkdir -p /opt/references && \
    wget -q https://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz -O /opt/references/hg19ToHg38.over.chain.gz && \
    wget -q https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes -O /opt/references/hg38.chrom.sizes

# Set environment variables
ENV LIFTOVER_CHAIN_PATH=/opt/references/hg19ToHg38.over.chain.gz
ENV CHROM_SIZES_PATH=/opt/references/hg38.chrom.sizes

LABEL maintainer="joanneboysen@gmail.com"