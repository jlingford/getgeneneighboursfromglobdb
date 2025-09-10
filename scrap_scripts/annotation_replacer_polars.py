#!/usr/bin/env python3

import polars as pl


# file paths
annotation_table = "../data/AMXMAG_0088_annotations_decompressed.tsv"
gff_input = "../data/AMXMAG_0088___3888___genetic_neighbourhood.gff"
new_gff_output = "../data/annotation_replacement_test3.gff"

# load master annotation table of genome
df = pl.read_csv(
    annotation_table,
    separator="\t",
    skip_lines=1,
    has_header=False,
    new_columns=["ID", "db_xref", "Name", "product", "evalue"],
    schema_overrides={"evalue": pl.Float64},
)

# keep only rows of KOfam annotations, and only the best evalue per gene ID
df = df.filter(pl.col("db_xref") == "KOfam")
df = df.sort("evalue").group_by("ID").first()


# load .gff file of genome
gff = pl.read_csv(
    gff_input,
    separator="\t",
    has_header=False,
    new_columns=[
        "gene_id",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ],
)

# gene ID is conatined within the attributes column. extract it, and use it as the key for merging with the annotations table
gff = gff.with_columns(pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID"))

# merge gff table and annotations table based on the shared gene ID
merged_gff = gff.join(df, on="ID", how="left")


# remake the gff attributes column using KOfam annotation information from the annotations table
merged_gff = merged_gff.with_columns(
    pl.when(pl.col("Name").is_not_null())
    .then(
        pl.format(
            "ID={};Name={};db_xref={};product={}",
            pl.col("ID"),
            pl.col("Name"),
            pl.col("db_xref"),
            pl.col("product"),
        )
    )
    .otherwise(pl.col("attributes"))
    .alias("attributes")
)

# discard the columns of the merged gff df that aren't native to the .gff format
new_gff = merged_gff.select(
    [
        "gene_id",
        "source",
        "type",
        "start",
        "end",
        "score",
        "strand",
        "phase",
        "attributes",
    ]
)

# write new gff file to output
new_gff.write_csv(
    new_gff_output,
    separator="\t",
    include_header=False,
)
