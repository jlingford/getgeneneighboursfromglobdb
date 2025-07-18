#!/usr/bin/env python3
"""
Extract and plot genetic neighbourhood around gene of interest
"""

# imports
import os
import re
import sys
import gzip
import shutil
import logging
import argparse
import textwrap
import numpy as np
import pandas as pd
import seaborn as sns
from Bio import SeqIO
from pathlib import Path
import matplotlib.pyplot as plt
import dna_features_viewer as dfv


def parse_arguments():
    """Parse arguments to script"""
    parser = argparse.ArgumentParser(
        description=textwrap.dedent("""\
            # -----------------------------#
            # Gene neighbourhood extractor #
            #------------------------------#
            A python script for finding gene neighbourhoods surrounding a target gene of interest from .gff files, and plotting them.
            Can extract multiple targets in bulk provided a list of target genes in a plain text file.
            --------------------------------
        """),
        epilog=textwrap.dedent("""\
            Examples:
            1) find gene neighbourhoods around a list of target genes, 5000 bp upstream, and 20,000 bp downstream of target genes. Plus add annotations:

            python %(prog)s -d path/to/gff_files_dir -l target_genes.txt -o output_dir -U 5000 -D 20000 -a path/to/annotation_dir
        """),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Add mutually exclusive groups of arguments
    mutparse = parser.add_mutually_exclusive_group(required=True)

    mutparse.add_argument(
        "-l",
        "--genelist",
        dest="gene_list",
        type=Path,
        metavar="FILE",
        help="Path to list of genes of interest (Either -l or -n are required).",
    )

    mutparse.add_argument(
        "-n",
        "--name",
        dest="gene_name",
        type=str,
        metavar="STR",
        help="Name of gene of interest (Either -l or -n are required)",
    )

    # Add individual arguments
    parser.add_argument(
        "-d",
        "--gffdir",
        dest="gff_path",
        type=Path,
        required=True,
        metavar="FILE|DIR",
        help="Path to .gff file OR path to base of dir containing .gff files (required)",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        dest="output_dir",
        type=str,
        default=".",
        required=False,
        metavar="DIR",
        help="Path to output dir. [Default: current dir]",
    )

    parser.add_argument(
        "-U",
        "--upstream",
        dest="upstream_window",
        type=int,
        default=10000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) upstream of GOI [Default: 10000]",
    )

    parser.add_argument(
        "-D",
        "--downstream",
        dest="downstream_window",
        type=int,
        default=10000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) downstream of GOI [Default: 10000]",
    )

    parser.add_argument(
        "-a",
        "--annodir",
        dest="annotation_dir",
        type=Path,
        required=False,
        metavar="DIR",
        help="Path to dir containing annotation files",
    )

    parser.add_argument(
        "-f",
        "--fastadir",
        dest="fasta_dir",
        type=Path,
        required=False,
        metavar="DIR",
        help="Path to dir containing fasta files",
    )

    # Parse arguments into args object
    args = parser.parse_args()

    # Validate arguments
    # if not args.gff_path.is_file():
    #     parser.error(f"Input file does not exist or is not a file: {args.input_dir}")

    # Info about input files
    if args.gff_path.is_file():
        print(f"Input .gff file provided: {args.gff_path}")

    if args.gff_path.is_dir():
        print(f"Target dir for .gff files provided. Searching {args.gff_path}/...")

    return args


def extract_gene_neighbourhood(args, gene_name, gff_file):
    """Finds the genes upstream and downstream of gene of interest in .gff file and returns/outputs new subset of the .gff file"""
    # params
    gff_file = gff_file  # The input gff file
    output_dir = args.output_dir  # Your output file
    gene_name = gene_name
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

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

    # Find gene of interest in .gff file
    # NOTE: hardcoded to look for ID=
    gene_row = gff[gff["attributes"].str.contains(f"ID={gene_name}", na=False)]

    if gene_row.empty:
        raise ValueError(f"Gene '{gene_name}' not found in GFF.")

    # Get genetic neighbourhood coordinates to center on gene of interest
    gene_start = gene_row.iloc[0]["start"]
    gene_end = gene_row.iloc[0]["end"]
    scaffold = gene_row.iloc[0]["seqid"]

    # Define window coordinates
    window_start = max(gene_start - upstream_window, 0)
    window_end = gene_end + downstream_window

    # Get subset of gff file based on window
    subset = gff[
        (gff["seqid"] == scaffold)
        & (gff["start"] <= window_end)
        & (gff["end"] >= window_start)
    ]

    # Write output
    outpath = Path(output_dir) / f"{gene_name}.gff"
    if outpath.exists():
        logging.warning(f"Writing over existing output .gff file: {outpath.name}")
    outpath.parent.mkdir(exist_ok=True, parents=True)
    subset.to_csv(
        outpath,
        sep="\t",
        header=False,
        index=False,
    )

    print(
        f"Extracted region ({window_start}-{window_end}) on {scaffold} written to: {outpath.name}"
    )


def plot_gene_neighbourhood():
    """Stuff"""


def fasta_neighbourhood_extract():
    """Extract fasta seqs surrounding gene of interest and reformat for AlphaPulldown input"""


def find_target_file(input_path: Path, target_name: str) -> Path:
    """Find target .gff file in input directory based on name of gene of interest"""
    matched_file = list(input_path.rglob(f"{target_name}*.gff*"))

    if not matched_file:
        raise FileNotFoundError(f"Could not find target file with name: {target_name}.")

    elif len(matched_file) > 1:
        logging.warning(f"Multiple files named: {target_name}: {matched_file}")
        logging.warning(f"Using first matched .gff file: {matched_file[0]}")

    return matched_file[0]


def process_target_genes(args):
    """Stuff"""
    # Passing list of genes of interest with specific .gff target file
    if args.gene_list and args.gff_path.is_file():
        with open(args.gene_list) as target_file:
            for line in target_file:
                gene_name = line.rstrip()
                # print(f"Finding gene neighbourhood around target gene: {gene_name}")
                extract_gene_neighbourhood(args, gene_name, args.gff_path)

    # Passing list of genes of interest with general gff parent dir
    if args.gene_list and args.gff_path.is_dir():
        with open(args.gene_list) as target_file:
            for line in target_file:
                gene_name = line.rstrip()
                target_name = gene_name.split("___")[0]
                gff_file = find_target_file(args.gff_path, target_name)
                extract_gene_neighbourhood(args, gene_name, gff_file)

    # Passing singular gene of interest with specific .gff target file
    if args.gene_name and args.gff_path.is_file():
        gene_name = args.gene_name
        extract_gene_neighbourhood(args, gene_name, args.gff_path)

    # Passing singular genes of interest with general gff parent dir
    if args.gene_name and args.gff_path.is_dir():
        gene_name = args.gene_name
        target_name = gene_name.split("___")[0]
        gff_file = find_target_file(args.gff_path, target_name)
        extract_gene_neighbourhood(args, gene_name, gff_file)


def main():
    args = parse_arguments()
    process_target_genes(args)


if __name__ == "__main__":
    main()
