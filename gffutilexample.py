#!/usr/bin/env python3

import pandas as pd

# --- Config ---
gff_path = "./AMXMAG_0088_cog_decompressed.gff"  # Your input GFF file
output_path = "nife_test.gff"  # Your output file
# gene_name = "AMXMAG_0088___3611"  # Replace with your gene of interest
gene_name = "AMXMAG_0088___3673"
window = 10000  # Size of window upstream/downstream

# --- Load GFF file into DataFrame ---
gff = pd.read_csv(
    data_dir,
    sep="\t",
    comment="#",
    header=None,
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

# --- Find the gene of interest by Name=... ---
# You can adjust the attribute search here if needed (e.g., ID=...)
goi_row = gff[gff["attributes"].str.contains(f"ID={gene_name}", na=False)]

if goi_row.empty:
    raise ValueError(f"Gene '{gene_name}' not found in GFF.")

# --- Get region coordinates ---
gene_start = goi_row.iloc[0]["start"]
gene_end = goi_row.iloc[0]["end"]
chrom = goi_row.iloc[0]["seqid"]

window_start = max(gene_start - window, 0)
window_end = gene_end + window

# --- Subset GFF for overlapping features on the same chromosome ---
subset = gff[
    (gff["seqid"] == chrom)
    & (gff["start"] <= window_end)
    & (gff["end"] >= window_start)
]

# --- Write to output ---
subset.to_csv(output_path, sep="\t", header=False, index=False)

print(
    f"Extracted region ({window_start}-{window_end}) on {chrom} written to: {output_path}"
)
