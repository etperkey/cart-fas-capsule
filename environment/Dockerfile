# FAS/FASLG CAR-T bystander killing — reproducibility environment.
# Pins to the versions reported in Methods: DESeq2 1.38.3, limma 3.54.2,
# GSVA 1.46.0, tximport 1.26.1, survival 3.4.0, survminer 0.4.9. Those
# Bioconductor releases correspond to Bioc 3.16 / R 4.2.3.
#
# The capsule reproduces bulk-RNA-seq Fig 1 and Sup Fig 1A–F only. The
# Haradhvala scRNA thread that Cooper's Box project explored was dropped
# from the final manuscript, so Seurat 5 / BPCells / SeuratDisk are NOT
# installed here.

FROM rocker/r-ver:4.2.3

ENV DEBIAN_FRONTEND=noninteractive \
    R_HOME=/usr/local/lib/R \
    TZ=America/Chicago

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libfontconfig1-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libfreetype6-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libcairo2-dev \
        libxt-dev \
        libglpk-dev \
        libgmp-dev \
        libgsl-dev \
        libbz2-dev \
        liblzma-dev \
        zlib1g-dev \
        pandoc \
    && rm -rf /var/lib/apt/lists/*

# Pin CRAN to a Posit Package Manager snapshot that corresponds to Bioc 3.16.
RUN R -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-04-20'), download.file.method = 'libcurl'); \
         install.packages(c('BiocManager','remotes','renv','rmarkdown','knitr','tidyverse','vroom','readxl','ggsci','ggpubr','ggtext','circlize','survival','survminer','maxstat','rstatix','patchwork','ggrepel','cowplot','scales','matrixStats','Matrix','R.utils'), Ncpus=4)"

RUN R -e "BiocManager::install(version = '3.16', ask = FALSE, update = FALSE); \
         BiocManager::install(c('DESeq2','edgeR','limma','tximport','GSVA','clusterProfiler','fgsea','GEOquery','AnnotationDbi','org.Hs.eg.db','GenomicFeatures','SummarizedExperiment','rtracklayer','ComplexHeatmap'), update = FALSE, ask = FALSE, Ncpus=4)"

# Extras pulled in by Cooper's `plotting_func.R` / `analysis_func.R` beyond
# the manuscript-declared Methods list. ggfortify + ggnewscale from CRAN
# (Posit snapshot pin); vsn from Bioc 3.16.
RUN R -e "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/jammy/2023-04-20'), download.file.method = 'libcurl'); \
         install.packages(c('ggfortify','ggnewscale'), Ncpus=4); \
         BiocManager::install('vsn', update = FALSE, ask = FALSE, Ncpus=4)"

WORKDIR /code
