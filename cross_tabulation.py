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