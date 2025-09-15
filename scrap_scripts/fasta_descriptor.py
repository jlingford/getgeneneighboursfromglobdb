#!/usr/bin/env python3

import gzip
from pathlib import Path
import polars as pl
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from io import TextIOWrapper


def open_fasta(file_path: Path) -> TextIOWrapper:
    """Open a fasta file that might be gzipped or plain text"""
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")


output_dir = "../data"
fasta_file = Path("../data/GCA_020628865.faa")
annotation_file = "../data/GCA_020628865_annotations.tsv"
# gene_ids = gene_ids_full
# gene_ids_nomaturation = gene_ids_nomaturation
# gene_ids_ppi_candidates = gene_ids_ppi_candidates
# gene_name = gene_name

# map annotation descriptions to gene ids
df = pl.read_csv(annotation_file, separator="\t", has_header=True, quote_char=None)
print(df)
df = df.select(pl.col("gene_callers_id"), pl.col("function"))
print(df)
df = df.unique(subset=["gene_callers_id"], keep="first")
print(df)
desc_map = df.transpose(column_names="gene_callers_id").to_dict(as_series=False)

print(desc_map)
# print(desc_map)
print(desc_map.get("GCA_020628865___1942", [None])[0])

# output all fastas from the entire gene neighbourhood
outpath_all = Path(output_dir) / "descriptive_fastas.faa"
outpath_all.parent.mkdir(exist_ok=True, parents=True)
with open_fasta(fasta_file) as handle:
    records = list(SeqIO.parse(handle, "fasta"))
    # fasta_subset = [rec for rec in records if rec.id in gene_ids]
    descriptive_fastas = [
        SeqRecord(
            seq=rec.seq,
            name=rec.id,
            id=rec.id,
            description=str(desc_map.get(rec.id, [None])[0]),
            # description="",
        )
        for rec in records
    ]
    with open(outpath_all, "w") as f:
        SeqIO.write(descriptive_fastas, f, "fasta")
