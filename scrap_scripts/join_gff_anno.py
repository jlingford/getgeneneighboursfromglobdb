#!/usr/bin/env python3

from operator import ge
import polars as pl
import re

gff_file = "../data/GCA_020628865_cog.gff"
annotation_file = "../data/GCA_020628865_annotations.tsv"
# gene_name = "GCA_020628865___632"
# gene_name = "GCA_020628865___125"
gene_name = "GCA_020628865___14"

# load master annotation table of genome
df = pl.read_csv(
    annotation_file,
    separator="\t",
    quote_char=None,
    has_header=True,
    new_columns=["ID", "db_xref", "Name", "product", "evalue"],
    schema_overrides={"evalue": pl.Float64},
)

print(df)

# keep only rows of KOfam or Pfam annotations, and only the best evalue per gene ID
# if args.use_kofam_annotation is True:
#     df = df.filter(pl.col("db_xref") == "KOfam")
# if args.use_pfam_annotation is True:
#     df = df.filter(pl.col("db_xref") == "Pfam")

# df = df.sort("evalue").group_by("ID").first()
# df = df.unique(subset=["ID"], maintain_order=False)

# print(df)

# load .gff file of genome
gff = pl.read_csv(
    gff_file,
    separator="\t",
    quote_char=None,
    comment_prefix="#",
    has_header=True,
    new_columns=[
        "seqid",
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

print(gff)

# gene ID is contained within the attributes column. extract it, and use it as the key for merging with the annotations table
gff = gff.with_columns(pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID"))

print(gff)

# merge gff table and annotations table based on the shared gene ID
merged_gff = gff.join(df, on="ID", how="left")

print(merged_gff)
# print(merged_gff.head(10).glimpse())

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

print(merged_gff)

# discard the columns of the merged gff df that aren't native to the .gff format
new_gff = merged_gff.select(
    [
        "seqid",
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

print(new_gff)
# print(new_gff.select(pl.col("attributes")).glimpse())

# outfile = "../data/test/new_gff.tsv"
# new_gff.write_csv(outfile, separator="\t")


pattern = rf"ID={re.escape(gene_name)}(?:;|$)"
# goi_row = gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}"))
goi_row = new_gff.filter(pl.col("attributes").str.contains(pattern))

print(goi_row)

# Get genetic neighbourhood coordinates to center on gene of interest
goi_start = goi_row[0, "start"]
goi_end = goi_row[0, "end"]
goi_scaffold = goi_row[0, "seqid"]
goi_strand = goi_row[0, "strand"]

# Define window coordinates
upstream_window = 10000
downstream_window = 10000
window_start: int = max(goi_start - upstream_window, 0)
window_end: int = goi_end + downstream_window

# Get subset of gff file based on window. Get genes from both strands, unless -S|--one_strand flag is provided
# gff of all genes in the target neighbourhood
gff_subset = new_gff.filter(
    (pl.col("strand") == goi_strand)
    & (pl.col("seqid") == goi_scaffold)
    & (pl.col("start") <= window_end)
    & (pl.col("end") >= window_start)
)

print(gff_subset)

SUBSET_LIST = [
    "K00635",  # 3
    "COG4908",  # 3
    "COG0456",  # 10
    "COG2267",  # 12
    "K03574",  # 13
    "COG1051",  # 13
]

gff_subset_from_list = new_gff.filter(
    (pl.col("strand") == goi_strand)
    & (pl.col("seqid") == goi_scaffold)
    & (pl.col("start") <= window_end)
    & (pl.col("end") >= window_start)
    & (pl.col("attributes").str.contains_any(SUBSET_LIST))
)

print(gff_subset_from_list)


outfile = "../data/test/new_gff_subset.tsv"
gff_subset.write_csv(outfile, separator="\t")
