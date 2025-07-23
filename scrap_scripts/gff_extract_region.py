#!/usr/bin/env python3

import pandas as pd

# --- Config ---
gff_path = "./AMXMAG_0088_cog_decompressed.gff"
output_path = "gff_test.gff"
# goi_name = "AMXMAG_0088___3611"
goi_name = "AMXMAG_0088___3673"
window_upstream = 10000
window_downstream = 10000

# --- Load GFF ---
gff = pd.read_csv(
    gff_path,
    sep="\t",
    comment="#",
    header=None,
    names=[
        "scaffold",
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

# Find gene of interest (GOI)
goi_row = gff[gff["attributes"].str.contains(f"ID={goi_name}", na=False)]

if goi_row.empty:
    raise ValueError(f"Gene '{goi_name}' not found in GFF.")

# Gene of interest position in scaffold
goi_start = goi_row.iloc[0]["start"]
goi_end = goi_row.iloc[0]["end"]
scaffold = goi_row.iloc[0]["scaffold"]

# Find scaffold length
scaffold_len = gff["end"].max()
print(scaffold)

# Set window coordinates in gff
window_start = goi_start - window_upstream
window_end = goi_end + window_downstream
# print(window_start, window_end)
total_window_size = window_end - window_start

# Get subset of gff file, only selecting for scaffold with GOI and with neighbouring genes that are within the window start and end
gff_subset = gff[
    (gff["scaffold"] == scaffold)
    & (gff["end"] >= window_start)
    & (gff["start"] <= window_end)
]

print(gff_subset)

# --- Write output ---
# gff_subset.to_csv(output_path, sep="\t", header=False, index=False)

print(
    f"Extracted genetic neighbourhood around GOI: {goi_name}. Region: {window_start} - {window_end} bp (total window size: {total_window_size}). Scaffold: {scaffold} "
)
