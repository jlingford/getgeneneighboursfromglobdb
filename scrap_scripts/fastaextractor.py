import pandas as pd
import gzip
from Bio import SeqIO

#######
# creating gff_subset

gff_path = "./data/AMXMAG_0088_cog.gff.gz"
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
chrom = goi_row.iloc[0]["seqid"]

window_start = max(gene_start - window, 0)
window_end = gene_end + window

# --- Subset GFF for overlapping features on the same chromosome ---
gff_subset = gff[
    (gff["seqid"] == chrom)
    & (gff["start"] <= window_end)
    & (gff["end"] >= window_start)
]

# extract gene_ids from .gff subset
gene_ids = gff_subset["attributes"].str.findall(r"ID=([^;]+)")
# create list of gene_ids from gene neighbourhood, getting them out of the nested list
gene_ids = [item for sublist in gene_ids for item in sublist]
print(gene_ids)

###################################################################

# extract fasta seqs from fasta file

fasta_file = "./data/AMXMAG_0088.faa.gz"
# fasta_file = "./data/AMXMAG_0088.faa"
gene_ids = gene_ids

with gzip.open(fasta_file, "rt") as handle:
    all_fastas = SeqIO.parse(
        handle,
        "fasta",
    )
    gene_neighbours = [record for record in all_fastas if record.id in gene_ids]
    with open("./output/subset_fasta.faa", "w") as f:
        SeqIO.write(gene_neighbours, f, "fasta")


# gene_neighbours = [record.seq for record in SeqIO.parse(file, "fasta") if record.name]
# fasta_subset = SeqIO.to_dict(gene_neighbours)

# for record in gene_neighbours:
#     print(">" + record.id + "\n" + record.seq)


# gene_neighbours = [
#     fasta_subset[gene_id] for gene_id in gene_ids if gene_id in fasta_subset
# ]
# print(gene_neighbours)


########################################################################
