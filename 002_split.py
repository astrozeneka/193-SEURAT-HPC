
import json
from glob import glob
import pandas as pd
import shapely
from shapely.geometry import shape
from os.path import basename

meta_csv = "larc_datasets/larc_merged_meta_data.csv"

if __name__ == '__main__':
    meta_df = pd.read_csv(meta_csv, index_col=0)
    x = meta_df["x_slide_mm"].values
    y = meta_df["y_slide_mm"].values

    for geojson_file in glob("geojson/*.geojson"):
        print(f"Processing {geojson_file}...")
        with open(geojson_file) as f:
            polygon = shape(json.load(f)["geometry"])
        slug = basename(geojson_file).replace(".geojson", "")

        mask = shapely.contains_xy(polygon, x, y)
        inside = meta_df[mask]
        inside.to_csv(f"splitted/meta/{slug}.csv")
        print(f"{geojson_file}: {mask.sum()} / {len(meta_df)} points inside")

    print("Done")