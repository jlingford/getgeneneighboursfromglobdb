#!/usr/bin/env python3

from operator import ge
from numpy import extract
from pandas import concat, merge
import polars as pl
import re

gff_file = "../data/GCA_020628865_cog.gff"
annotation_file = "../data/GCA_020628865_annotations.tsv"
# gene_name = "GCA_020628865___632"
gene_name = "GCA_020628865___125"
# gene_name = "GCA_020628865___14"
# gene_name = "GCA_020628865___41"

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

# print(gff)

# gene ID is contained within the attributes column. extract it, and use it as the key for merging with the annotations table
gff2 = gff.with_columns(pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID"))

# print(gff)

# merge gff table and annotations table based on the shared gene ID
merged_df = gff2.join(df, on="ID", how="left")

# print(merged_df)
# print(merged_df.head(10).glimpse())

# remake the gff attributes column using KOfam annotation information from the annotations table
merged_gff = merged_df.with_columns(
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
).select(
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

# print("new merged gff:")
# print(merged_gff)

# # discard the columns of the merged gff df that aren't native to the .gff format
# new_gff = merged_gff.select(
#     [
#         "seqid",
#         "source",
#         "type",
#         "start",
#         "end",
#         "score",
#         "strand",
#         "phase",
#         "attributes",
#     ]
# )

# print(new_gff)
# print(new_gff.select(pl.col("attributes")).glimpse())

# outfile = "../data/test/new_gff.tsv"
# new_gff.write_csv(outfile, separator="\t")


# pattern = rf"ID={re.escape(gene_name)}(?:;|$)"
# pattern = f"ID={gene_name}(?:;|$)"
# goi_row = gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}"))
# goi_row = merged_gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}(?:;|$)"))
#
# print("goi row:")
# print(goi_row)
#
# # Get genetic neighbourhood coordinates to center on gene of interest
# goi_start = goi_row[0, "start"]
# goi_end = goi_row[0, "end"]
# goi_scaffold = goi_row[0, "seqid"]
# goi_strand = goi_row[0, "strand"]
#
# # Define window coordinates
# upstream_window = 10000
# downstream_window = 10000
# window_start: int = max(goi_start - upstream_window, 0)
# window_end: int = goi_end + downstream_window
#
# # Get subset of gff file based on window. Get genes from both strands, unless -S|--one_strand flag is provided
# # gff of all genes in the target neighbourhood
# gff_subset = merged_gff.filter(
#     (pl.col("strand") == goi_strand)
#     & (pl.col("seqid") == goi_scaffold)
#     & (pl.col("start") <= window_end)
#     & (pl.col("end") >= window_start)
# )
#
# print(gff_subset)
#
# SUBSET_LIST = [
#     "K00635",  # 3
#     "COG4908",  # 3
#     "COG0456",  # 10
#     "COG2267",  # 12
#     "K03574",  # 13
#     "COG1051",  # 13
# ]
#
# gff_subset_from_list = merged_gff.filter(
#     (pl.col("strand") == goi_strand)
#     & (pl.col("seqid") == goi_scaffold)
#     & (pl.col("start") <= window_end)
#     & (pl.col("end") >= window_start)
#     & (pl.col("attributes").str.contains_any(SUBSET_LIST))
# )
#
# print(gff_subset_from_list)
#
#
# gff_subset_from_list: list = (
#     gff_subset_from_list.select(
#         pl.col("attributes")
#         .str.extract_groups(r"ID=([^;]+)")
#         .struct.field("1")
#         .alias("gene_id")
#     )
#     .to_series()
#     .unique()
#     .to_list()
# )
#
# print(gff_subset_from_list)
#
# # outfile = "../data/test/new_gff_subset.tsv"
# # gff_subset.write_csv(outfile, separator="\t")
#
#
# # Figuring out how to filter merged_gff by db_xref:
#
# print(merged_df)
#
# # gf = merged_df.with_columns(
# #     pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("geneID"),
# #     pl.col("attributes").str.extract(r"db_xref=([^;]+)", 1).alias("HMM_type"),
# # )

gf = merged_df.select(
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
        "ID",
        "db_xref",
        "evalue",
    ]
)

# print(gf)


gf = gf.with_columns(
    pl.col("ID").str.extract(r"___(\d+)$", 1).cast(pl.Int32).alias("ID_num")
)

gf = gf.sort(by=["ID_num", "evalue"])

print("sorted gf")
print(gf)

# uniq_gf = gf.unique(subset=["ID_num"], keep="first").sort(by=["ID_num"])
#
# print("unique gf:")
# print(uniq_gf)

# filt = (
#     gf.filter((pl.col("db_xref").is_null()) | (pl.col("db_xref") == "Pfam"))
#     .group_by("ID_num", maintain_order=True)
#     .first()
# ).sort(by="ID_num")
#
# print(filt)
#
# filt2 = (
#     gf.filter((pl.col("db_xref").is_null()) | (pl.col("db_xref") == "COG20_FUNCTION"))
#     .group_by("ID_num", maintain_order=True)
#     .first()
# ).sort(by="ID_num")
#
# print(filt2)
#
# filt3 = (
#     gf.filter((pl.col("db_xref").is_null()) | (pl.col("db_xref") == "KOfam"))
#     .group_by("ID", maintain_order=True)
#     .first()
# ).sort(by="ID_num")
#
# print(filt3)
# filt.write_csv("../data/test/new_gff_subset2.tsv", separator="\t")


# gff without any info after first attributes, and only the ID in attributes
clean_gff = gff.with_columns(pl.col("attributes").str.extract(r"(ID=[^;]+)", 1))

print("clean_gff:")
print(clean_gff)

# gff with new columns of ID (join key) and ID_num (for sorting)
clean_gff1 = gff.with_columns(
    pl.col("attributes").str.extract(r"(ID=[^;]+)", 1),
    pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID"),
    pl.col("attributes")
    .str.extract(r"___(\d+)(?:;|$)", 1)
    .cast(pl.Int32)
    .alias("ID_num"),
)

print("clean_gff1:")
print(clean_gff1)

# clean_gff2 = clean_gff.with_columns(
#     pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID")
# )
# print(clean_gff2)

merged_df = clean_gff1.join(df, on="ID", how="left")

# print("merged_df:")
# print(merged_df)

# merge1: make keep db_xref to clean_gff
merged_gff1 = merged_df.sort(by=["ID_num", "evalue"]).select(
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
        # "db_xref",
        "ID_num",
    ]
)

print("clean merged_gff1:")
print(merged_gff1)

# merge2: populate gff attributes with annotation info
merged_gff2 = (
    merged_df.with_columns(
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
    # .select(
    #     [
    #         "seqid",
    #         "source",
    #         "type",
    #         "start",
    #         "end",
    #         "score",
    #         "strand",
    #         "phase",
    #         "attributes",
    #         "ID",
    #         "db_xref",
    #         "evalue",
    #     ]
    # )
    # .with_columns(
    #     pl.col("ID").str.extract(r"___(\d+)$", 1).cast(pl.Int32).alias("ID_num")
    # )
    .sort(by=["ID_num", "evalue"])
    .select(
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
            # "db_xref",
            "ID_num",
        ]
    )
)

print("merged_gff2:")
print(merged_gff2)

concat_gf = pl.concat([merged_gff2, merged_gff1]).sort(by="ID_num")

print("concat_gf:")
print(concat_gf)

final_gf = (
    concat_gf.unique(subset="attributes").sort(by="ID_num").select(pl.exclude("ID_num"))
)

# print("final_gf:")
# print(final_gf)
#
# expand_gf = final_gf.with_columns(
#     pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID"),
#     pl.col("attributes").str.contains(r"db_xref=([^;]+)").alias("any_dbxref"),
#     pl.col("attributes").str.contains(r"db_xref=Pfam").alias("has_pfam"),
#     pl.col("attributes").str.contains(r"db_xref=KOfam").alias("has_kofam"),
#     pl.col("attributes")
#     .str.contains(r"db_xref=COG20_FUNCTION")
#     .alias("has_cog20function"),
# )
#
# print("expand_gf:")
# print(expand_gf)


# filt = (
#     final_gf.filter((pl.col("db_xref").is_null()) | (pl.col("db_xref") == "KOfam"))
#     .group_by("ID_num", maintain_order=True)
#     .first()
# ).sort(by="ID_num")
#
# print(filt)
