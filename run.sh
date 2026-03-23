#!/bin/bash
# run container named 193-seurat-hpc
docker run -it \
    -v ./larc_datasets:/app/larc_datasets \
    -v ./cell_list:/app/cell_list \
    -v ./stats:/app/stats \
    193-seurat-hpc:latest
