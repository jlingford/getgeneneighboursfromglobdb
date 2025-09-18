#!/usr/bin/env python3

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


# load .gff file of genome
raw_gff = pl.read_csv(
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

pfam_annodf = full_annodf.filter(pl.col("db_xref") == "Pfam")
kofam_annodf = full_annodf.filter(pl.col("db_xref") == "KOfam")
cogfun_annodf = full_annodf.filter(pl.col("db_xref") == "COG20_FUNCTION")

gff = raw_gff.with_columns(
    pl.col("attributes").str.extract(r"ID=([^;]+)(?:;|$)", 1).alias("ID"),
    pl.col("attributes")
    .str.extract(r"ID=.*___(\d+)(?:;|$)", 1)
    .cast(pl.Int32)
    .alias("ID_num"),
)

# print(gff)

full_merge = gff.join(full_annodf, on="ID", how="left").sort(by=["ID_num", "evalue"])
pfam_merge = gff.join(pfam_annodf, on="ID", how="left").sort(by=["ID_num", "evalue"])
kofam_merge = gff.join(kofam_annodf, on="ID", how="left").sort(by=["ID_num", "evalue"])
cogfun_merge = gff.join(cogfun_annodf, on="ID", how="left").sort(
    by=["ID_num", "evalue"]
)

# print(full_merge)

full_popgff = (
    full_merge.with_columns(
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
            "evalue",
            "ID",
            "ID_num",
        ]
    )
    .sort(by=["ID_num", "evalue"])
    .unique(subset="ID", keep="first")
    .sort(by="ID_num")
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
        ]
    )
)

print(full_popgff)


pfam_popgff = (
    pfam_merge.with_columns(
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
            "evalue",
            "ID",
            "ID_num",
        ]
    )
    .sort(by=["ID_num", "evalue"])
    .unique(subset="ID", keep="first")
    .sort(by="ID_num")
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
        ]
    )
)

print(pfam_popgff)
