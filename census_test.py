import cellxgene_census

with cellxgene_census.open_soma(census_version="2025-11-08") as census:

    # Pull that dataset's human cell metadata into a pandas dataframe
    df = cellxgene_census.get_obs(
        census,
        "homo_sapiens",
        value_filter="tissue_general == 'amygdala'",
        column_names=["cell_type", "tissue", "tissue_general"],
    )


print("Shape:", df.shape)
print(df["cell_type"].value_counts().head(10))