#!/usr/bin/env python3

from numpy import printoptions
import pandas as pd

gff_file = "../data/AMXMAG_0088___3888___genetic_neighbourhood.gff"
gene_name = "AMXMAG_0088___3887"
upstream_window = 5000
downstream_window = 5000

# --- Load GFF file into DataFrame ---
gff = pd.read_csv(
    gff_file,
    sep="\t",
    comment="#",
    header=None,
    compression="infer",
    names=[
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

# Find gene of interest (GOI) in .gff file
# NOTE: hardcoded to look for ID=
goi_row = gff[gff["attributes"].str.contains(f"ID={gene_name}", na=False)]
if goi_row.empty:
    raise ValueError(f"Gene '{gene_name}' not found in target .gff file.")

# Get genetic neighbourhood coordinates to center on gene of interest
goi_start = goi_row.iloc[0]["start"]
goi_end = goi_row.iloc[0]["end"]
goi_scaffold = goi_row.iloc[0]["seqid"]
goi_strand = goi_row.iloc[0]["strand"]

# Define window coordinates
window_start: int = max(goi_start - upstream_window, 0)
window_end: int = goi_end + downstream_window

# Get subset of gff file based on window. Get only genes from same strand as GOI, unless specified otherwise
gff_subset = gff[
    (gff["strand"] == goi_strand)
    & (gff["seqid"] == goi_scaffold)
    & (gff["start"] <= window_end)
    & (gff["end"] >= window_start)
]

#     gff_subset = gff[
#         (gff["seqid"] == goi_scaffold)
#         & (gff["start"] <= window_end)
#         & (gff["end"] >= window_start)
#     ]
# get subset of gff without maturation genes, for smaller fasta file extraction
# TODO: incorporate this function in a smarter way:
gff_subset_nomaturation = gff[
    (gff["strand"] == goi_strand)
    & (gff["seqid"] == goi_scaffold)
    & (gff["start"] <= window_end)
    & (gff["end"] >= window_start)
    & (~gff["attributes"].str.contains("maturation", na=False))
]

# # Write output
# outpath = Path(output_dir) / gene_name / f"{gene_name}___genetic_neighbourhood.gff"
# if outpath.exists():
#     logging.warning(f"Writing over existing output .gff file: {outpath.name}")
# outpath.parent.mkdir(exist_ok=True, parents=True)
# gff_subset.to_csv(
#     outpath,
#     sep="\t",
#     header=False,
#     index=False,
# )
# print(
#     f"Extracted region ({window_start}-{window_end}) on {goi_scaffold} written to: {outpath.name}"
# )

# match gene ID name contained within "ID="
gene_ids = gff_subset["attributes"].str.findall(r"ID=([^;]+)")
# create list of gene_ids from gene neighbourhood, getting them out of the nested list
gene_ids: list = [item for sublist in gene_ids for item in sublist]

# get subset of gff without maturation genes, for smaller fasta file extraction
# TODO: incorporate this function in a smarter way:
gene_ids_nomaturation = gff_subset_nomaturation["attributes"].str.findall(r"ID=([^;]+)")
gene_ids_nomaturation: list = [
    item for sublist in gene_ids_nomaturation for item in sublist
]

print(gene_ids)
print(gene_ids_nomaturation)
