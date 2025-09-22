#!/usr/bin/env python3

import polars as pl

# get hydrogenase classification
file = "../../nife_geneids2hyddb.tsv"
gene_name = "MOTU40_070268___637"
# gene_name = "GCF_026891035___2818"
# gene_name = "AMXMAG_0088___3888"

hyddb_df = pl.read_csv(
    file,
    separator="\t",
    has_header=False,
    new_columns=["gene_id", "classification"],
)

print(hyddb_df)

# hydclass = hyddb_df.filter(
#     pl.when(pl.col("gene_id") == f"{gene_name}").select("classification").item()
# )

hydclass = (
    hyddb_df.filter(pl.col("gene_id") == f"{gene_name}")
    .get_column("classification")
    .first()
)


# hydclass = hyddb_df.filter(
#     pl.when(pl.col("gene_id") == f"{gene_name}").is_not_null().then(
#         pl.select("classification").item()
#     )
# )

print(hydclass)
