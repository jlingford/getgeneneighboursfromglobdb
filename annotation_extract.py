import pandas as pd
from pandas.core import missing


# --- Config ---
gff_path = "./data/AMXMAG_0088_cog_decompressed.gff"  # Your input GFF file
output_path = "nife_test.gff"  # Your output file
gene_name = "AMXMAG_0088___3611"  # Replace with your gene of interest
# gene_name = "AMXMAG_0088___3673"
window = 10000  # Size of window upstream/downstream

# --- Load GFF file into DataFrame ---
gff = pd.read_csv(
    gff_path,
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
scaffold = goi_row.iloc[0]["seqid"]

window_start = max(gene_start - window, 0)
window_end = gene_end + window

# --- Subset GFF for overlapping features on the same chromosome ---
gff_subset = gff[
    (gff["seqid"] == scaffold)
    & (gff["start"] <= window_end)
    & (gff["end"] >= window_start)
]

#############################################################
# params
gff_subset = gff_subset
anno_file = "./data/AMXMAG_0088_annotations_decompressed.tsv"
gene_of_interest = "AMXMAG_0088___3611"
anno = pd.read_csv(anno_file, delimiter="\t", compression="infer")

# Reformat annotation table
# match gene ID name contained within "ID="
gene_ids = gff_subset["attributes"].str.findall(r"ID=([^;]+)")

# create list of gene_ids, getting them out of the nested list
gene_ids = [item for sublist in gene_ids for item in sublist]
# print(gene_ids)

# locate row containing gene of interest
goi_row = anno[anno["gene_callers_id"].str.contains(gene_of_interest)]

# create subset of annotation table based on gene_ids list, and sorted based on gene_ids and evalues
anno_subset = anno[anno["gene_callers_id"].isin(gene_ids)].sort_values(
    by=["gene_callers_id", "e_value"], ascending=[True, True]
)

# find anny missing annotations based on gene_ids
anno_ids = anno_subset["gene_callers_id"].unique()
missing_ids = set(gene_ids) - set(anno_ids)
print(missing_ids)

gene_ids_df = pd.DataFrame({"gene_callers_id": gene_ids})
# print(gene_ids_df)

# Merge: keep all gene_ids, add data from final_result where available
merged_with_missing = gene_ids_df.merge(anno_subset, on="gene_callers_id", how="left")

# Optionally fill missing values for empty rows
merged_filled = merged_with_missing.fillna(
    {
        "source": float("nan"),
        "accession": float("nan"),
        "function": float("nan"),
        "e_value": float("nan"),
    }
)
print(merged_filled)

merged_filled.to_csv(outpath, sep="\t", index=False)


##################################################
# getting subset of just KOfam annotations and next best annotations... more complicated

# get table of kofam annotations, but lowest evalue per unique gene_id
kofam_rows = anno_subset[anno_subset["source"].str.contains("KOfam")].loc[
    lambda df: df.groupby("gene_callers_id")["e_value"].idxmin()
]

# Identify gene_callers_ids that have no KOfam entry
all_ids = anno_subset["gene_callers_id"].unique()
kofam_ids = kofam_rows["gene_callers_id"].unique()
missing_ids = set(all_ids) - set(kofam_ids)
# print(missing_ids)

# For those missing, take best overall row (lowest e_value)
fallback_rows = anno_subset[anno_subset["gene_callers_id"].isin(missing_ids)].loc[
    lambda df: df.groupby("gene_callers_id")["e_value"].idxmin()
]

# Combine both
merged_df = (
    pd.concat([kofam_rows, fallback_rows])
    .sort_values("gene_callers_id")
    .reset_index(drop=True)
)

# print(final_result)

gene_ids_df = pd.DataFrame({"gene_callers_id": gene_ids})
# print(gene_ids_df)

# Merge: keep all gene_ids, add data from final_result where available
merged_with_missing = gene_ids_df.merge(merged_df, on="gene_callers_id", how="left")

# Optionally fill missing values for empty rows
merged_filled = merged_with_missing.fillna(
    {
        "source": float("nan"),
        "accession": float("nan"),
        "function": float("nan"),
        "e_value": float("nan"),  # Or use "-" if preferred
    }
)

# print(merged_filled)


##############################################################################
# anno_subset.to_csv("./output/anno_subset.tsv", sep="\t", index=False)

# df_melted = anno_subset.melt(
#     id_vars=["gene_callers_id", "source"],
#     value_vars=["accession", "function", "e_value"],
#     var_name="field",
#     value_name="value",
# )
#
# print(df_melted)

# Ensure e_value is numeric
# anno_subset["e_value"] = pd.to_numeric(anno_subset["e_value"], errors="coerce")
#
# gene_name["gene_callers_id"].sort_values()

# get table of kofam annotations
# kofam_rows = anno_subset.loc[anno_subset["source"].str.contains("KOfam")]
# print(kofam_rows)
# kofam_rows = anno_subset.loc[anno_subset["source"].str.contains("KOfam")]
# print(kofam_rows)
# kofam_rows = anno_subset[anno_subset["source"] == "KOfam"]
# print(kofam_rows)
#
# merged_anno = gene_ids_df.merge(anno_subset, on="gene_callers_id", how="left")
# print(merged_anno)
# non_kofam_rows = anno_subset[anno_subset["source"] != "KOfam"]
#
# get table of lowerst evalue rows
# best_rows = anno_subset.loc[anno_subset.groupby("gene_callers_id")["e_value"].idxmin()]
