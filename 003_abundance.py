from glob import glob
import pandas as pd
from os.path import basename


if __name__ == '__main__':
    relative_abundance_stacked = None # aka percentage
    for file in glob("splitted/meta/*.csv"):
        print(f"Processing {file}...")
        df = pd.read_csv(file, index_col=0)
        df_cell_types = df[["RNA_Standard_Cell.Typing.InSituType.1_1_clusters"]]
        sample_name = basename(file).split("_")[0]
        # Compute percentage in a stacked way
        cell_type_counts = df_cell_types.value_counts().sort_index()
        cell_type_percentage = cell_type_counts / cell_type_counts.sum()
        cell_type_percentage.name = sample_name
        if relative_abundance_stacked is None:
            relative_abundance_stacked = cell_type_percentage
        else:
            relative_abundance_stacked = pd.concat([relative_abundance_stacked, cell_type_percentage], axis=1)
    # Save the relative abundance to csv
    relative_abundance_stacked.to_csv("splitted/overall_abundance/overall_abundance.csv")
    print("done")