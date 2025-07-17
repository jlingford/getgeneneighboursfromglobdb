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
        description="Script for gene neighbourhoods surrounding a target gene of interest from .gff files, and plotting them. Can extract multiple targets in bulk.",
        epilog="Example: python %(prog)s -d <path/to/gff_files_dir> -t <path/to/target_gene_list.txt> [-o output_dir] [-U INT] [-D INT] [-a annotation_dir]",
    )

    # Add arguments
    parser.add_argument(
        "-t",
        "--target-gene-list",
        dest="target_gene_list",
        type=Path,
        required=True,
        help="Path to input file (required).",
    )

    parser.add_argument(
        "-d",
        "--input-dir",
        dest="input_dir",
        type=Path,
        required=True,
        help="Path to base of dir containing .gff files (required)",
    )

    parser.add_argument(
        "-o",
        "--output-dir",
        dest="output_dir",
        type=Path,
        default=Path("."),
        required=False,
        help="Path to output dir. [Default: current dir]",
    )

    parser.add_argument(
        "-U",
        "--upstream-window",
        dest="upstream_window",
        type=int,
        default=10000,
        required=False,
        help="Size of window (in base pairs) upstream of GOI [Default: 10000]",
    )

    parser.add_argument(
        "-D",
        "--downstream-window",
        dest="downstream_window",
        type=int,
        default=10000,
        required=False,
        help="Size of window (in base pairs) downstream of GOI [Default: 10000]",
    )

    parser.add_argument(
        "-a",
        "--annotation-dir",
        dest="annotation_dir",
        type=Path,
        required=False,
        help="Path to dir containing annotation files",
    )

    # Parse arguments into args object
    args = parser.parse_args()

    # Validate arguments
    if not args.input_file.is_file():
        parser.error(f"Input file does not exist or is not a file: {args.input_file}")

    if not args.output_file.is_dir():
        parser.error(
            f"Output directory does not exist or is not a directory: {args.output_dir}"
        )

    return args


def extract_gene_neighbourhood():
    """Stuff"""


def plot_gene_neighbourhood():
    """Stuff"""


def process_target_genes():
    """Stuff"""


def func(args):
    """Stuff goes here"""
    # do stuff...


def main():
    args = parse_arguments()
    primary_function(args)


if __name__ == "__main__":
    main()
