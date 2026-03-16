FROM rocker/r-ver:4.1.1

# System dependencies for Seurat and HDF5
RUN apt-get update && apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN Rscript -e "\
    install.packages(c('ggplot2', 'dplyr'), repos='https://cloud.r-project.org'); \
    install.packages('BiocManager', repos='https://cloud.r-project.org'); \
    BiocManager::install(version='3.14'); \
    install.packages('Seurat', repos='https://cloud.r-project.org') \
"

CMD ["R"]
