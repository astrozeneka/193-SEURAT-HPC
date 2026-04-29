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
RUN Rscript -e "install.packages('Seurat', repos='https://cloud.r-project.org');"
RUN Rscript -e "BiocManager::install('dorothea'); \
    BiocManager::install('viper'); "

RUN Rscript -e "install.packages('devtools'); \
    install.packages('remotes')";


RUN Rscript -e "library('remotes'); \
    devtools::install_github(repo = 'hhoeflin/hdf5r'); \
    devtools::install_github(repo = 'mojaveazure/loomR', ref = 'develop')"
RUN Rscript -e "remotes::install_github( \
      'Nanostring-Biostats/CosMx-Analysis-Scratch-Space', \
      subdir = '_code/scPearsonPCA', \
      ref = 'Main' \
    ); \
    install.packages(c('data.table', 'ggplot2', 'ggrepel'));  \
    install.packages('pals'); \
"

RUN Rscript -e "install.packages('harmony');"

RUN Rscript -e "library('remotes'); \
    install.packages('Matrix', type = 'source'); \
    install.packages('irlba', type = 'source'); \
    devtools::install_github('https://github.com/Nanostring-Biostats/InSituType'); \
"

RUN apt-get update && apt-get install -y libglpk-dev && rm -rf /var/lib/apt/lists/*


RUN Rscript -e "library('remotes'); \
    print('Hello!'); \
    devtools::install_github('saeyslab/nichenetr')\
"

RUN Rscript -e "install.packages('leidenbase')"


RUN Rscript -e "BiocManager::install(c(\
  'BiocNeighbors',\
  'ComplexHeatmap',\
  'BiocParallel'\
));\
install.packages(c(\
  'igraph', 'ggplot2', 'ggalluvial', 'svglite',\
  'Seurat', 'reticulate', 'umap', 'ggnetwork',\
  'FNN', 'sna', 'network', 'pbapply', 'future',\
  'future.apply', 'RColorBrewer', 'NMF', 'circlize',\
  'ggpubr', 'patchwork', 'Matrix'\
))"

RUN apt-get update && apt-get install -y libuv1-dev && rm -rf /var/lib/apt/lists/*


RUN Rscript -e "library('remotes'); \
    devtools::install_github('jinworks/CellChat')"

WORKDIR /app

CMD ["R"]
