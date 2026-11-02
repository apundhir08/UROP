import cellxgene_census

cols = [
    "cell_type",
    "dataset_id",
    "disease",
    "tissue",
    "tissue_general",
    "assay",
    "suspension_type",
    "is_primary_data",
]

with cellxgene_census.open_soma(census_version="2025-11-08") as census:
    df = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter="tissue == 'large intestine'",
        column_names=cols,
    )

print("Before dedup:", df.shape)
print(df["is_primary_data"].value_counts())

df_all = df.copy()  
df = df[df["is_primary_data"] == True]


print("After dedup:", df.shape)

print("Number of datasets:", df["dataset_id"].nunique())
counts = df["dataset_id"].value_counts()
print(counts[counts > 0])

print("\nCell types:")
print(df["cell_type"].value_counts().head(10))

d = df["disease"].value_counts()
print(d[d > 0])

#Print de-duplicated version
print("\nPrimary data:")
print(df["is_primary_data"].value_counts())

import pandas as pd

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", None)

def report(data, label):
    data = data.copy()
    for c in ["cell_type", "disease", "dataset_id"]:
        data[c] = data[c].astype(str)   # avoids zero-count categories
    short = data["dataset_id"].str[:8]

    print(f"===== {label}: {data.shape[0]} cells, {data['dataset_id'].nunique()} datasets =====")

    top = data["cell_type"].value_counts().head(5).index
    sub = data[data["cell_type"].isin(top)]
    print("Cell type x dataset:")
    print(pd.crosstab(sub["cell_type"], sub["dataset_id"].str[:8]), "\n")

    print("Disease x dataset (cell counts):")
    print(pd.crosstab(data["disease"], short), "\n")

    print("Disease x dataset (fraction of each dataset's cells):")
    print(pd.crosstab(short, data["disease"], normalize="index").round(2), "\n")

report(df_all, "All datasets (before dedup)")
report(df, "Primary data only")