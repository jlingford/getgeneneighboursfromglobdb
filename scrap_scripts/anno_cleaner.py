#!/usr/bin/env python3

from numpy import full
import polars as pl


gff_file = "../data/GCA_020628865_cog.gff"
annotation_file = "../data/GCA_020628865_annotations.tsv"

# load master annotation table of genome
full_annodf = pl.read_csv(
    annotation_file,
    separator="\t",
    quote_char=None,
    has_header=True,
    new_columns=["ID", "db_xref", "Name", "product", "evalue"],
    schema_overrides={"evalue": pl.Float64},
)


# print(full_annodf)

# ha = full_annodf.filter(pl.col("Name").str.contains("!!!"))
# print(ha)

clean_anno = full_annodf.with_columns(
    pl.col("Name").str.extract(r"([^!]+)"), pl.col("product").str.extract(r"([^!]+)")
)

print(clean_anno)

ha = clean_anno.filter(pl.col("ID").str.contains(r"GCA_020628865___7$"))
print(ha)
