from pathlib import Path
from Bio import SeqIO
import polars as pl


target_gene = "GCA_020628865___632"

codes_of_interest = [
    "COG1740",  # Ni,Fe-hydrogenase I small subunit (HyaA) (PDB:6FPI)
    "COG1035",  # Coenzyme F420-reducing hydrogenase, beta subunit (FrhB) (PDB:3ZFS)
    "PF14720.10",  # NiFe/NiFeSe hydrogenase small subunit C-terminal
    "K06441",  # ferredoxin hydrogenase gamma subunit [EC:1.12.7.2]
    "K06282",  # hydrogenase small subunit [EC:1.12.99.6]
    "K00441",  # coenzyme F420 hydrogenase subunit beta [EC:1.12.98.1]
    "K14113",  # energy-converting hydrogenase B subunit D
    "K14127",  # F420-non-reducing hydrogenase iron-sulfur subunit [EC:1.12.99.- 1.8.98.5 1.8.98.6]
    "K18006",  # [NiFe] hydrogenase diaphorase moiety small subunit [EC:1.12.1.2]
    "K17992",  # NADP-reducing hydrogenase subunit HndB [EC:1.12.1.3]
    "K18006",  # [NiFe] hydrogenase diaphorase moiety small subunit [EC:1.12.1.2]
    "K23548",  # uptake hydrogenase small subunit [EC:1.12.99.6]
    "K05927",  # quinone-reactive Ni/Fe-hydrogenase small subunit [EC:1.12.5.1]
]

gff_file = "../data/GCA_020628865_cog.gff"
fasta_file = "../data/GCA_020628865.faa"
output_fasta = "../data/test/SSUfasta.faa"
output_fasta2 = "../data/test/LSU-SSU-fasta.faa"
upstream_window = 2000
downstream_window = 2000

# --- Load GFF ---
gff = pl.read_csv(
    gff_file,
    separator="\t",
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

# --- Locate target gene ---
goi = gff.filter(pl.col("attributes").str.contains(f"ID={target_gene}"))

scaffold = goi.item(0, "seqid")
strand = goi.item(0, "strand")
start = goi.item(0, "start")
end = goi.item(0, "end")

window_start = max(start - upstream_window, 0)
window_end = end + downstream_window

# --- Extract neighbourhood ---
neighbours = gff.filter(
    (pl.col("seqid") == scaffold)
    & (pl.col("start") <= window_end)
    & (pl.col("end") >= window_start)
    & (pl.col("strand") == strand)
)

# --- Identify candidate matching codes ---
candidates = neighbours.filter(
    pl.any_horizontal(
        [pl.col("attributes").str.contains(code) for code in codes_of_interest]
    )
)

if candidates.is_empty():
    raise ValueError(
        f"No neighbours near {target_gene} match codes {codes_of_interest}"
    )

# Choose the closest neighbour (e.g., smallest distance to target)
candidates = candidates.with_columns(
    pl.when(pl.col("start") > end)
    .then(pl.col("start") - end)
    .otherwise(start - pl.col("end"))
    .alias("distance")
)
closest = candidates.sort("distance").head(1)

# --- Extract gene ID from attributes ---
neighbour_id = (
    closest.select(
        pl.col("attributes").str.extract_groups(r"ID=([^;]+)").struct.field("1")
    )
    .to_series()
    .to_list()[0]
)

print(neighbour_id)
print(target_gene)

# --- Filter FASTA for neighbour ---
records = list(SeqIO.parse(fasta_file, "fasta"))
matching_record = next((rec for rec in records if rec.id == neighbour_id), None)
target_record = next((rec for rec in records if rec.id == target_gene), None)
combined_records = [target_record, matching_record]

print("matching_record:")
print(matching_record)
print("")
print("target_record:")
print(target_record)
print(combined_records)

# --- Write to FASTA ---
SeqIO.write(matching_record, output_fasta, "fasta")
print(f"Neighbour {neighbour_id} saved to {output_fasta}")

SeqIO.write(combined_records, output_fasta2, "fasta")
print(f"{target_gene} saved to {output_fasta}")
