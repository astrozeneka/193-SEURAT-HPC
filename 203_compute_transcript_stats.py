import pandas as pd

if __name__ == '__main__':

    counts = pd.read_csv("stats/transcript_counts.csv")

    stats = pd.DataFrame({
        "metric": ["median_n_transcripts", "median_n_unique_transcripts"],
        "value": [
            counts["n_transcripts"].median(),
            counts["n_unique_transcripts"].median()
        ]
    })

    stats.to_csv("stats/transcript_stats.csv", index=False)
    print("Done")
