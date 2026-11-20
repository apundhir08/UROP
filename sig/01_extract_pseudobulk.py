import sys; sys.path.insert(0, ".")
import numpy as np, pandas as pd, scipy.sparse as sp
import cellxgene_census
from census_conn import tiledb_cfg

VERSION   = "2025-11-08"
DATASET   = "d18736c3-6292-4379-919a-d6d973204c87"
CELL_TYPE = "classical monocyte"
OBS_COLS  = ["donor_id", "disease", "sex", "development_stage"]

with cellxgene_census.open_soma(census_version=VERSION, tiledb_config=tiledb_cfg()) as census:
    adata = cellxgene_census.get_anndata(
        census, organism="Homo sapiens", X_name="raw",
        obs_value_filter=(f"dataset_id == '{DATASET}' and cell_type == '{CELL_TYPE}' "
                          "and is_primary_data == True"),
        obs_column_names=OBS_COLS)

X = adata.X.tocsr()
assert np.allclose(X.data, np.round(X.data)), "raw layer is not integer counts"
obs = adata.obs.assign(donor_id=adata.obs["donor_id"].astype(str))

codes, donors = pd.factorize(obs["donor_id"])
agg = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                    shape=(len(donors), len(codes)))
pb = (agg @ X).toarray().astype(np.int64)

donor_tbl = (obs.groupby("donor_id", observed=True)
                .agg(disease=("disease", "first"), sex=("sex", "first"),
                     stage=("development_stage", "first"),
                     n_cells=("donor_id", "size"), n_labels=("disease", "nunique"))
                .reindex(donors))
assert (donor_tbl["n_labels"] == 1).all() and donor_tbl["n_cells"].sum() == adata.n_obs
donor_tbl["lib_size"] = pb.sum(axis=1)

genes = adata.var[["feature_id", "feature_name"]].reset_index(drop=True)
pd.DataFrame(pb.T, index=genes["feature_id"], columns=donors).to_csv("sig/pseudobulk_counts.csv")
genes.to_csv("sig/genes.csv", index=False)
donor_tbl.to_csv("sig/donor_table.csv")
print("cells:", adata.n_obs, "| genes:", adata.n_vars, "| donors:", len(donors))
print(donor_tbl["disease"].astype(str).value_counts().to_dict())