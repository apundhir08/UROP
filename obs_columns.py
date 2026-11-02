import cellxgene_census

with cellxgene_census.open_soma(census_version="2025-11-08") as census:

    # Pull that dataset's human cell metadata into a pandas dataframe
    obs = census["census_data"]["homo_sapiens"].obs
    columns = list(obs.keys())
print("Number of obs columns:", len(columns))
for i, name in enumerate(columns, start=1):
    print(f"{i:2d}. {name}")