import scipy.io
import pandas as pd
import numpy as np

# load MEX format sparse matrix
matrix = scipy.io.mmread("rna_mex/matrix.mtx").tocsr()
genes = pd.read_csv("rna_mex/genes.csv")["x"].tolist()
barcodes = pd.read_csv("rna_mex/barcodes.csv")["x"].tolist()

if __name__ == '__main__':

    # Total transcript counts per cell (sum of all UMI counts per column)
    counts_per_cell = np.asarray(matrix.sum(axis=0)).flatten()

    # Number of unique genes detected per cell (non-zero entries per column)
    matrix_csc = matrix.tocsc()
    unique_genes_per_cell = np.diff(matrix_csc.indptr)

    summary = pd.DataFrame({
        "cell_id": barcodes,
        "n_transcripts": counts_per_cell,
        "n_unique_transcripts": unique_genes_per_cell
    })

    summary.to_csv("stats/transcript_counts.csv", index=False)
    print("Done")