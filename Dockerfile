FROM rocker/r-ver:4.5

# System dependencies for Seurat and HDF5
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libpng-dev \
    libfftw3-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN Rscript -e "\
    install.packages(c('ggplot2', 'dplyr'), repos='https://cloud.r-project.org'); \
    install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install(version='3.22'); \
    BiocManager::install('multtest'); \
    install.packages('Seurat', repos='https://cloud.r-project.org'); \
    BiocManager::install('OmnipathR'); \
    BiocManager::install('decoupleR') \
"

RUN Rscript -e "BiocManager::install('DESeq2')"

RUN Rscipt -e "BiocManager::install('dorothea'); \
    BiocManager::install('viper'); "

RUN Rscript -e "remotes::install_github( \
      'Nanostring-Biostats/CosMx-Analysis-Scratch-Space', \
      subdir = '_code/scPearsonPCA', \
      ref = 'Main' \
    ); \
    install.packages(c('data.table', 'ggplot2', 'ggrepel'));  \
    install.packages('pals'); \
"

RUN Rscript -e "install.packages('harmony');"

WORKDIR /app

CMD ["R"]
