

# To fix docker api compatibility issues, set the following environment variable:
export DOCKER_API_VERSION=1.44


Merge objects (with slide ID)
    ↓
QC filter per slide
    ↓
Normalize (log-normalization)
    ↓
Assign group labels in metadata
    ↓
Extract normalized matrix → split by group
    ↓
Per gene: Log2FC + Wilcoxon p-value
    ↓
BH correction → adjusted p-value
    ↓
Volcano plot (Log2FC × -log10 adj.p)


ssh ryanr@node1