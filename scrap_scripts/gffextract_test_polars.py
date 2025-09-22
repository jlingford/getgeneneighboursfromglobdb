#!/usr/bin/env python3

from numpy import int64
import polars as pl
from pathlib import Path

gfffile = "/run/media/james/T7/globdb_r226/globdb_r226_gff_cog/R226CHUNK000/AMXMAG_0088_cog.gff.gz"
outfile = "./gff_test.tsv"
# gene_name = "AMXMAG_0088___3888"
gene_name = "AMXMAG_0088___3878"

# params
upstream_window = 100000
downstream_window = 100000

# step 1: get gff subset
raw_gff = pl.read_csv(
    gfffile,
    separator="\t",
    quote_char=None,
    comment_prefix="#",
    has_header=False,
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

print(raw_gff)

# Find gene of interest (GOI) in .gff file
goi_row = raw_gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}(?:;|$)"))

print(goi_row)

# # Get genetic neighbourhood coordinates to center on gene of interest
# # OLD WAY
# goi_start = goi_row[0, "start"]
# goi_end = goi_row[0, "end"]
# goi_scaffold = goi_row[0, "seqid"]
# goi_strand = goi_row[0, "strand"]

# Ensure scalar extraction
goi_start = goi_row.select("start").item()
goi_end = goi_row.select("end").item()
goi_scaffold = goi_row.select("seqid").item()
goi_strand = goi_row.select("strand").item()

goi_size = goi_end - goi_start

# print(f"goi start: {type(goi_start)}")
print(f"goi start: {goi_start}")
print(f"goi end: {goi_end}")
print(f"goi size: {goi_size}")
print(f"scaffold name: {goi_scaffold}")
print(f"goi strand: {goi_strand}")

# Define window coordinates, conditional on strand direction
if goi_strand == "+":
    window_start: int = max(goi_start - upstream_window, 1)
    window_end: int = goi_end + downstream_window
else:
    window_start: int = goi_start + upstream_window
    window_end: int = max(goi_end - downstream_window, 1)

# print(f"Window start: {type(window_start)}")
print(f"Window start: {window_start}")
print(f"Window end: {window_end}")

# Get subset of gff file based on window. Get genes from both strands, unless -S|--one_strand flag is provided
# if args.one_strand is True:
#     # gff of all genes in the target neighbourhood
#     gff_subset = raw_gff.filter(
#         (pl.col("strand") == goi_strand)
#         & (pl.col("seqid") == goi_scaffold)
#         & (pl.col("start") <= window_end)
#         & (pl.col("end") >= window_start)
#         & (pl.col("type") == "CDS")
#     )
# else:
# gff of all genes in the target neighbourhood
if goi_strand == "+":
    gff_subset = raw_gff.filter(
        (pl.col("seqid") == goi_scaffold)
        & (pl.col("start") <= window_end)
        & (pl.col("end") >= window_start)
        & (pl.col("type") == "CDS")
    )
else:
    gff_subset = raw_gff.filter(
        (pl.col("seqid") == goi_scaffold)
        & (pl.col("start") >= window_end)
        & (pl.col("end") <= window_start)
        & (pl.col("type") == "CDS")
    )


print(gff_subset)

# return list of gene IDs for fasta file output, called at end of this function
gene_ids: list = (
    gff_subset.select(
        pl.col("attributes")
        .str.extract_groups(r"ID=([^;]+)")
        .struct.field("1")
        .alias("gene_id")
    )
    .to_series()
    .unique()
    .to_list()
)

# remove other attributes info and add key column to merge and sort on
blank_gff = gff_subset.with_columns(
    pl.col("attributes").str.extract(r"ID=([^;]+)(?:;|$)", 1).alias("ID"),
    pl.col("attributes")
    .str.extract(r"ID=.*___(\d+)(?:;|$)", 1)
    .cast(pl.Int32)
    .alias("ID_num"),
)
