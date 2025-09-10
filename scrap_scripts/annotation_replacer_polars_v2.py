#!/usr/bin/env python3

import polars as pl
import re


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

gene_name = "AMXMAG_0088___3887"

# NOTE: hardcoded to look for ID=
goi_row = gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}"))
# if goi_row.empty:
#     raise ValueError(f"Gene '{gene_name}' not found in target .gff file.")

# Get genetic neighbourhood coordinates to center on gene of interest
goi_start = goi_row[0, "start"]
goi_end = goi_row[0, "end"]
goi_scaffold = goi_row[0, "gene_id"]
goi_strand = goi_row[0, "strand"]

window_size = 5000

# Define window coordinates
window_start: int = max(goi_start - window_size, 0)
window_end: int = goi_end + window_size

# Get subset of gff file based on window. Get only genes from same strand as GOI, unless specified otherwise
gff_subset = gff.filter(
    (pl.col("strand") == goi_strand)
    & (pl.col("gene_id") == goi_scaffold)
    & (pl.col("start") <= window_end)
    & (pl.col("end") >= window_start)
)
gff_subset = gff.filter(
    (pl.col("gene_id") == goi_scaffold)
    & (pl.col("start") <= window_end)
    & (pl.col("end") >= window_start)
)
# get subset of gff without maturation genes, for smaller fasta file extraction
# TODO: incorporate this function in a smarter way:
gff_subset_nomaturation = gff.filter(
    (pl.col("strand") == goi_strand)
    & (pl.col("gene_id") == goi_scaffold)
    & (pl.col("start") <= window_end)
    & (pl.col("end") >= window_start)
    & (~pl.col("attributes").str.contains("maturation"))
)

gene_ids: list = (
    gff_subset.select(
        pl.col("attributes")
        .str.extract_groups(r"ID=([^;]+)")
        .struct.field("1")
        .alias("gene_id")
    )
    .to_series()
    .to_list()
)

gene_ids_nomaturation: list = (
    gff_subset_nomaturation.select(
        pl.col("attributes")
        .str.extract_groups(r"ID=([^;]+)")
        .struct.field("1")
        .alias("gene_id")
    )
    .to_series()
    .to_list()
)

test = (
    gff_subset.select(
        pl.col("attributes")
        .str.extract_groups(r"ID=([^;]+)")
        .struct.field("1")
        .alias("gene_id")
    )
    .to_series()
    .to_list()
)


print(gene_ids)
print(gene_ids_nomaturation)
