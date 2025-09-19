#!/usr/bin/env python3

from numpy import dtype, select
import polars as pl
from pathlib import Path
import pickle


def find_taxonomy_file_from_index(file_index: dict[str, list[Path]]) -> Path:
    """use prebuilt index to find GlobDB taxonomy file"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == r"globdb_r226_taxonomy.tsv"
        or filename == r"globdb_r226_taxonomy.tsv.gz"
        for path in paths
    ]
    if not matches:
        raise filenotfounderror("ERROR: No taxonomy file found")
    return matches[0]


index_path = Path("../globdb_index.pkl")

if index_path.is_file():
    with open(index_path, "rb") as f:
        pickle_index = pickle.load(f)

taxonomy_file = find_taxonomy_file_from_index(pickle_index)


def taxonomy_dataframe(taxonomy_file: Path) -> pl.DataFrame:
    """Read GlobDB taxonomy file and return dataframe of taxa"""
    # read in tsv file
    df = pl.read_csv(
        taxonomy_file,
        separator="\t",
        has_header=False,
        new_columns=["genome_id", "taxonomy_string"],
    )
    # remove prefix from taxa names, split taxonomy_string into new columns on ";", and rename columns to taxa
    df = (
        df.with_columns(
            pl.col("taxonomy_string")
            .str.replace_all(r"[a-z]__", "")
            .str.split_exact(by=";", n=6)
            .alias("taxa")
        )
        .unnest("taxa")
        .rename(
            {
                f"field_{i}": name
                for i, name in enumerate(
                    ["domain", "phylum", "class", "order", "family", "genus", "species"]
                )
            }
        )
    )
    return df


df = taxonomy_dataframe(taxonomy_file)
result = df.filter(pl.col("genome_id") == "MOTU40_070268").select("species").item()

print(str(result))
