
import pandas as pd
from glob import glob
from os.path import basename

CLINICAL_DATA_PATH = '../125-WHOLE-SLIDE-ANALYSIS/density_with_cd/Clinical_data.csv'
META_DIR = 'splitted/meta/*.csv'


if __name__ == '__main__':

    clinical_df = pd.read_csv(CLINICAL_DATA_PATH)
    clinical_df = clinical_df.dropna(subset=['SampleId'])
    clinical_df.set_index('SampleId', inplace=True)

    patient_response = clinical_df['pCR'].to_dict()

    cells_by_response = {
        "Responder": set(),
        "Non-responder": set()
    }
    for file in glob(META_DIR):
        patient_id = basename(file).replace(".csv", "")
        if patient_id not in patient_response:
            print(f"WARNING: no clinical data for {patient_id} — skipping.")
            continue
        response = patient_response[patient_id]
        cell_ids = pd.read_csv(file, usecols=[0], index_col=0).index
        cells_by_response[response].update(cell_ids)

    # save to 'Responder.txt" and "Non-responder.txt"
    for response, cell_ids in cells_by_response.items():
        with open(f"cell_list/{response}.txt", "w") as f:
            for cell_id in cell_ids:
                f.write(f"{cell_id}\n")