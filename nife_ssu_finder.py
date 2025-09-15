#!/usr/bin/env python3
"""
Extract and plot genetic neighbourhood around gene of interest

Control flow:
>main()
    ↳ parse_arguments()
    ↳ process_target_genes()
        ↳ extract_gene_neighbourhood()
            ↳ replace_gff_attributes()
            ↳ plot_gene_neighbourhood()
            ↳ fasta_neighbourhood_extract()
                ↳ fasta_pairwise_generation()
        ↳ annotation_extract()
"""
# TODO:
# change -n to nargs=+
# add logging output to argparse function

# imports
from ast import List
from io import TextIOWrapper
from os import path
from string import ascii_uppercase
import pickle
import gzip
import logging
import argparse
from operator import ge
import textwrap
from numpy import argsort
import polars as pl
import seaborn as sns
from BCBio import GFF
from pathlib import Path
from typing import TextIO
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
from itertools import combinations_with_replacement


# WARN: hardcoded HMM annotation codes for NiFe LSU
NIFE_LSU_CODES = [
    "COG0374",  # Ni,Fe-hydrogenase I large subunit (HyaB) (PDB:6FPI) (PUBMED:25905665)
    "COG0649",  # NADH:ubiquinone oxidoreductase 49 kD subunit (chain D) (NuoD) (PDB:2FUG)
    "COG3259",  # Coenzyme F420-reducing hydrogenase, alpha subunit (FrhA) (PDB:3ZFS)
    "COG3261",  # Ni,Fe-hydrogenase III large subunit (HycE2) (PDB:6CFW)
    "COG4042",  # Energy-converting hydrogenase Eha subunit A (EhaA) (PUBMED:19495416)
    "K00333",  # NADH-quinone oxidoreductase subunit D [EC:7.1.1.2]
    "K00338",  # NADH-quinone oxidoreductase subunit I [EC:7.1.1.2]
    "K00436",  # NAD-reducing hydrogenase large subunit [EC:1.12.1.2]
    "K00440",  # coenzyme F420 hydrogenase subunit alpha [EC:1.12.98.1]
    "K05588",  # bidirectional [NiFe] hydrogenase diaphorase subunit [EC:7.1.1.2]
    "K06281",  # hydrogenase large subunit [EC:1.12.99.6]
    "K13378",  # NADH-quinone oxidoreductase subunit C/D [EC:7.1.1.2]
    "K14086",  # ech hydrogenase subunit A
    "K14126",  # F420-non-reducing hydrogenase large subunit [EC:1.12.99.- 1.8.98.5]
    "K23549",  # uptake hydrogenase large subunit [EC:1.12.99.6]
    "PF00346.23",  # Respiratory-chain NADH dehydrogenase, 49 Kd subunit
    "PF00346.23",  # Respiratory-chain NADH dehydrogenase, 49 Kd subunit
    "PF00374.23",  # Nickel-dependent hydrogenase
    "PF17367.6",  # NiFe-hydrogenase-type-3 Eha complex subunit A
]

# WARN: hardcoded HMM annotation codes for NiFe SSU
NIFE_SSU_CODES = [
    "COG0377",  # NADH:ubiquinone oxidoreductase 20 kD subunit (chain B) or related Fe-S oxidoreductase (NuoB) (PDB:2FUG) (PUBMED:25941396)
    "COG1740",  # Ni,Fe-hydrogenase I small subunit (HyaA) (PDB:6FPI)
    "COG1941",  # Coenzyme F420-reducing hydrogenase, gamma subunit (FrhG) (PDB:5ODC)
    "COG3260",  # Ni,Fe-hydrogenase III small subunit (HycG) (PDB:6CFW)
    "K00331",  # NADH-quinone oxidoreductase subunit B [EC:7.1.1.2]
    "K05927",  # quinone-reactive Ni/Fe-hydrogenase small subunit [EC:1.12.5.1]
    "K06282",  # hydrogenase small subunit [EC:1.12.99.6]
    "K06441",  # ferredoxin hydrogenase gamma subunit [EC:1.12.7.2]
    "K14113",  # energy-converting hydrogenase B subunit D
    "K14127",  # F420-non-reducing hydrogenase iron-sulfur subunit [EC:1.12.99.- 1.8.98.5 1.8.98.6]
    "K14128",  # F420-non-reducing hydrogenase small subunit [EC:1.12.99.- 1.8.98.5]
    "K18006",  # [NiFe] hydrogenase diaphorase moiety small subunit [EC:1.12.1.2]
    "K18007",  # NAD-reducing hydrogenase small subunit [EC:1.12.1.2]
    "K23548",  # uptake hydrogenase small subunit [EC:1.12.99.6]
    "PF01058.26",  # NADH ubiquinone oxidoreductase, 20 Kd subunit
    "PF14720.10",  # NiFe/NiFeSe hydrogenase small subunit C-terminal
]

# WARN: hardcoded HMM annotation codes for NiFe FrhD / MvhD subunit (group 3 NiFe)
NIFE_FRHB_CODES = [
    "COG1035",  # Coenzyme F420-reducing hydrogenase, beta subunit (FrhB) (PDB:3ZFS)
    "COG1908",  # Coenzyme F420-reducing hydrogenase, delta subunit (MvhD) (PDB:5ODC)
    "K00441",  # coenzyme F420 hydrogenase subunit beta [EC:1.12.98.1]
    "K17992",  # NADP-reducing hydrogenase subunit HndB [EC:1.12.1.3]
    "PF02662.20",  # Methyl-viologen-reducing hydrogenase, delta subunit
]

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_GROUP1_PPI_CODES = [
    "COG1969",  # Ni,Fe-hydrogenase I cytochrome b subunit (HyaC) (PDB:4GD3)
    "COG3658",  # Cytochrome b subunit of Ni2+-dependent hydrogenase (CytB)
    "COG5557",  # Ni/Fe-hydrogenase 2 integral membrane subunit HybB (HybB)
    "K03620",  # Ni/Fe-hydrogenase 1 B-type cytochrome subunit
    "K04013",  # cytochrome c-type protein NrfB
    "PF01292.24",  # Prokaryotic cytochrome b561
    "PF14522.10",  # Cytochrome c7 and related cytochrome c
    "PF14537.10",  # Cytochrome c3
    "PF14537.10",  # Cytochrome c3
]

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_GROUP2_PPI_CODES = []

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_GROUP3_PPI_CODES = []

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_GROUP4_PPI_CODES = []

# WARN: hardcoded HMM annotation codes for all possible NiFe PPI partners
NIFE_PPI_CANDIDATES = (
    NIFE_LSU_CODES
    + NIFE_SSU_CODES
    + NIFE_FRHB_CODES
    + NIFE_GROUP1_PPI_CODES
    + NIFE_GROUP2_PPI_CODES
    + NIFE_GROUP3_PPI_CODES
    + NIFE_GROUP4_PPI_CODES
)

# WARN: hardcoded HMM annotation codes for NiFe group 1 partners
NIFE_MATURATION_CODES = [
    "COG0068",  # Hydrogenase maturation factor HypF (carbamoyltransferase) (HypF) (PDB:3TSP)
    "COG0298",  # Hydrogenase maturation factor HybG, HypC/HupF family (HypC) (PDB:2OT2)
    "COG0309",  # Carbamoyl dehydratase HypE (hydrogenase maturation factor) (HypE) (PDB:2I6R) (PUBMED:28618091)
    "COG0375",  # Hydrogenase maturation factor HypA/HybF, metallochaperone involved in Ni insertion (HybF) (PDB:3A43) (PUBMED:16817932)
    "COG0378",  # Hydrogenase/urease maturation factor HypB, Ni2+-binding GTPase (HypB) (PDB:2HF8) (PUBMED:24338018;23899293)
    "COG0409",  # Hydrogenase maturation factor HypD (HypD) (PDB:2Z1D)
    "COG0680",  # Ni,Fe-hydrogenase maturation factor (HyaD) (PDB:1CFZ) (PUBMED:15294295)
    "K03605",  # hydrogenase maturation protease [EC:3.4.23.-]
    "K03618",  # hydrogenase-1 operon protein HyaF
    "K04651",  # hydrogenase nickel incorporation protein HypA/HybF
    "K04652",  # hydrogenase nickel incorporation protein HypB
    "K04653",  # hydrogenase expression/formation protein HypC
    "K04654",  # hydrogenase expression/formation protein HypD
    "K04655",  # hydrogenase expression/formation protein HypE
    "K04656",  # hydrogenase maturation protein HypF
    "K07321",  # CO dehydrogenase maturation factor
    "PF01155.23",  # Hydrogenase/urease nickel incorporation, metallochaperone, hypA
    "PF01455.22",  # HupF/HypC family
    "PF01924.20",  # Hydrogenase formation hypA family
    "PF02492.23",  # CobW/HypB/UreG, nucleotide-binding domain
    "PF04809.17",  # HupH hydrogenase expression protein, C-terminal conserved region
    "PF11939.12",  # [NiFe]-hydrogenase assembly, chaperone, HybE
]


def parse_arguments() -> argparse.Namespace:
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

            python %(prog)s -d path/to/gff_files_dir -l target_genes.txt -o output_dir -U 5000 -D 20000 -a path/to/annotation_files
        """),
        formatter_class=argparse.RawTextHelpFormatter,
    )

    # Add mutually exclusive groups of arguments
    input_genes_args = parser.add_mutually_exclusive_group(required=True)
    globdb_input_args = parser.add_mutually_exclusive_group(required=True)
    anno_args = parser.add_mutually_exclusive_group()

    input_genes_args.add_argument(
        "-l",
        "--gene_list",
        dest="gene_list",
        type=Path,
        metavar="FILE",
        help="Path to list of genes of interest (Either -l or -n is required).",
    )

    input_genes_args.add_argument(
        "-n",
        "--gene_name",
        dest="gene_name",
        type=str,
        metavar="STR",
        help="Name of gene of interest (Either -l or -n is required)",
    )

    # Add individual arguments
    globdb_input_args.add_argument(
        "-d",
        "--globdb_dir",
        dest="data_dir",
        type=Path,
        metavar="DIR",
        help="Path to base of GlobDB dir containing all .gff, .tsv, and .faa files (Either -d or -i is required)",
    )

    globdb_input_args.add_argument(
        "-i",
        "--globdb_index",
        dest="index_path",
        type=Path,
        metavar="FILE.pkl",
        help="Path to precomputed .pkl index file of GlobDB directory structure (Either -d or -i is required)",
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
        default=5000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) upstream of GOI [Default: 5000]",
    )

    parser.add_argument(
        "-D",
        "--downstream",
        dest="downstream_window",
        type=int,
        default=5000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) downstream of GOI [Default: 5000]",
    )

    parser.add_argument(
        "--no_annotation",
        dest="no_annotation",
        action="store_true",
        help="Option to not return a table of annotations for the gene neighbourhood [Default: off]",
    )

    parser.add_argument(
        "--no_fasta",
        dest="no_fasta",
        action="store_true",
        help="Option to not return a fasta file for the gene neighbourhood [Default: off]",
    )

    parser.add_argument(
        "--no_ssu",
        dest="no_ssu",
        action="store_true",
        help="Option to not return a fasta file for the NiFe SSU [Default: off]",
    )

    parser.add_argument(
        "--no_plot",
        dest="no_plot",
        action="store_true",
        help="Option to not output genetic neighbourhood plots [Default: off]",
    )

    parser.add_argument(
        "--dpi",
        dest="plot_dpi",
        type=int,
        default=300,
        required=False,
        metavar="INT",
        help="Resolution of .png output plot in dpi [Default: 300]",
    )

    parser.add_argument(
        "--plot_format",
        dest="plot_format",
        choices=["png", "pdf", "svg"],
        default="png",
        required=False,
        metavar="STR",
        help="File format for output plot. Choices = png, pdf, svg [Default: png]",
    )

    parser.add_argument(
        "-B",
        "--both_strands",
        dest="both_strands",
        action="store_true",
        help="Option to include both strands of genetic neighbourhood in output, rather than only genetic neighbours on the same strand as the gene-of-interest. [Default: off]",
    )

    parser.add_argument(
        "--pairwise_set1",
        dest="pairwise_set1",
        action="store_true",
        help="Generate pairwise combinations of all sequences in gene neighbourhood, excluding maturation factors [Default: off]",
    )

    parser.add_argument(
        "--pairwise_set2",
        dest="pairwise_set2",
        action="store_true",
        help="Generate pairwise combinations of sequences in gene neighbourhood, but limited to most likely [NiFe] hydrogenase PPI candidates [Default: off]",
    )

    parser.add_argument(
        "--chai_fastas",
        dest="chai_fastas",
        action="store_true",
        help="Option to generate pairwise fasta combinations for Chai-1 input [Default: off]",
    )

    parser.add_argument(
        "--boltz_fastas",
        dest="boltz_fastas",
        action="store_true",
        help="Option to format pairwise fasta combinations for Boltz-1/2 input [Default: off]",
    )

    parser.add_argument(
        "--colabfold_fastas",
        dest="colabfold_fastas",
        action="store_true",
        help="Option to format pairwise fasta combinations for ColabFold search/batch input [Default: off]",
    )

    parser.add_argument(
        "--regular_fastas",
        dest="regular_fastas",
        action="store_true",
        help="Option add pairwise fasta combinations with plain fasta format [Default: on if --chai_fastas|boltz_fastas|colabfold_fastas options are not set]",
    )

    anno_args.add_argument(
        "-K",
        "--use_kofam",
        dest="use_kofam_annotation",
        action="store_true",
        help="Replace default COG20_FUNCTION annotations with KOfam annotations",
    )

    anno_args.add_argument(
        "-P",
        "--use_pfam",
        dest="use_pfam_annotation",
        action="store_true",
        help="Replace default COG20_FUNCTION annotations with Pfam annotations",
    )

    # Parse arguments into args object
    args = parser.parse_args()

    # check pairwise fasta args
    if (
        args.boltz_fastas
        or args.chai_fastas
        or args.colabfold_fastas
        or args.regular_fastas
    ) and not (args.pairwise_set1 or args.pairwise_set2):
        parser.error(
            "--pairwise_set1 or --pairwise_set2 options must be set too alongside the --boltz_fastas, --chai_fastas, --colabfold_fastas, and --regular_fastas options"
        )

    if (args.no_fasta) and (args.pairwise_set1 or args.pairwise_set2):
        parser.error("--no_fasta and --pairwise_set1/2 are mutually exclusive options")

    return args


def replace_gff_attributes(
    args: argparse.Namespace,
    gene_name: str,
    gff_file: Path,
    annotation_file: Path,
) -> pl.DataFrame:
    """Rewrite the attributes column of the gff file with KOfam annotation information instead of the default COG20_FUNCTION annotations"""
    # params
    gff_file = gff_file
    annotation_file = annotation_file
    gene_name = gene_name

    # load master annotation table of genome
    df = pl.read_csv(
        annotation_file,
        separator="\t",
        quote_char=None,
        has_header=True,
        new_columns=["ID", "db_xref", "Name", "product", "evalue"],
        schema_overrides={"evalue": pl.Float64},
    )

    # keep only rows of KOfam or Pfam annotations, and only the best evalue per gene ID
    if args.use_kofam_annotation is True:
        df = df.filter(pl.col("db_xref") == "KOfam")
    if args.use_pfam_annotation is True:
        df = df.filter(pl.col("db_xref") == "Pfam")

    # df = df.sort("evalue").group_by("ID").first()
    df = df.unique(subset=["ID"], maintain_order=True)

    # load .gff file of genome
    gff = pl.read_csv(
        gff_file,
        separator="\t",
        quote_char=None,
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

    # gene ID is conatined within the attributes column. extract it, and use it as the key for merging with the annotations table
    gff = gff.with_columns(
        pl.col("attributes").str.extract(r"ID=([^;]+)", 1).alias("ID")
    )

    # merge gff table and annotations table based on the shared gene ID
    merged_gff = gff.join(df, on="ID", how="left")

    # remake the gff attributes column using KOfam annotation information from the annotations table
    merged_gff = merged_gff.with_columns(
        pl.when(pl.col("Name").is_not_null())
        .then(
            pl.format(
                "ID={};Name={};db_xref={};product={}",
                pl.col("ID"),
                pl.col("Name"),
                pl.col("db_xref"),
                pl.col("product"),
            )
        )
        .otherwise(pl.col("attributes"))
        .alias("attributes")
    )

    # discard the columns of the merged gff df that aren't native to the .gff format
    new_gff = merged_gff.select(
        [
            "seqid",
            "source",
            "type",
            "start",
            "end",
            "score",
            "strand",
            "phase",
            "attributes",
        ]
    )

    return new_gff


def extract_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_file: Path,
    anno_file: Path,
    fasta_file: Path,
) -> None:
    """Finds the genes upstream and downstream of gene of interest in .gff file and returns/outputs new subset of the .gff file"""
    # params
    output_dir = args.output_dir  # Your output file
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

    # load gff dataframe
    if args.use_kofam_annotation is True:
        gff = replace_gff_attributes(args, gene_name, gff_file, anno_file)
    elif args.use_pfam_annotation is True:
        gff = replace_gff_attributes(args, gene_name, gff_file, anno_file)
    else:
        gff = pl.read_csv(
            gff_file,
            separator="\t",
            quote_char=None,
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

    # Find gene of interest (GOI) in .gff file
    # NOTE: hardcoded to look for ID=
    goi_row = gff.filter(pl.col("attributes").str.contains(f"ID={gene_name}"))
    # if goi_row.empty:
    #     raise ValueError(f"Gene '{gene_name}' not found in target .gff file.")

    # Get genetic neighbourhood coordinates to center on gene of interest
    goi_start = goi_row[0, "start"]
    goi_end = goi_row[0, "end"]
    goi_scaffold = goi_row[0, "seqid"]
    goi_strand = goi_row[0, "strand"]

    # Define window coordinates
    window_start: int = max(goi_start - upstream_window, 0)
    window_end: int = goi_end + downstream_window

    # Get subset of gff file based on window. Get only genes from same strand as GOI, unless -B|--both_strands flag is provided
    if args.both_strands is True:
        # gff of all genes in the target neighbourhood
        gff_subset = gff.filter(
            (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
        )
        # gff subset WITHOUT maturation genes, for smaller fasta file extraction
        gff_subset_nomaturation = gff.filter(
            (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
            & (~pl.col("attributes").str.contains_any(NIFE_MATURATION_CODES))
        )
        # gff subset WITH ONLY relevant PPI candidates for even smaller fasta file candidates
        gff_subset_ppi_candidates = gff.filter(
            (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
            & (pl.col("attributes").str.contains_any(NIFE_PPI_CANDIDATES))
        )
    else:
        # gff of all genes in the target neighbourhood
        gff_subset = gff.filter(
            (pl.col("strand") == goi_strand)
            & (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
        )
        # gff subset WITHOUT maturation genes, for smaller fasta file extraction
        gff_subset_nomaturation = gff.filter(
            (pl.col("strand") == goi_strand)
            & (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
            & (~pl.col("attributes").str.contains_any(NIFE_MATURATION_CODES))
        )
        # gff subset WITH ONLY relevant PPI candidates for even smaller fasta file candidates
        gff_subset_ppi_candidates = gff.filter(
            (pl.col("strand") == goi_strand)
            & (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
            & (~pl.col("attributes").str.contains_any(NIFE_PPI_CANDIDATES))
        )

    # return list of gene IDs for fasta file output, called at end of this function
    gene_ids: list = (
        gff_subset.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .to_list()
    )

    gene_ids_nomaturation: list = (
        gff_subset_nomaturation.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .to_list()
    )

    gene_ids_ppi_candidates: list = (
        gff_subset_ppi_candidates.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .to_list()
    )

    # Write output
    if args.use_kofam_annotation is True:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_KOfam_annotations.gff"
        )
    elif args.use_pfam_annotation is True:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_Pfam_annotations.gff"
        )
    else:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_COG20_annotations.gff"
        )

    if gff_outpath.exists():
        logging.warning(f"Writing over existing output .gff file: {gff_outpath.name}")
    gff_outpath.parent.mkdir(exist_ok=True, parents=True)

    gff_subset.write_csv(
        gff_outpath,
        separator="\t",
        include_header=False,
    )
    print(
        f"Extracted region ({window_start}-{window_end}) on {goi_scaffold} written to: {gff_outpath.name}"
    )

    # Execute downstream functions, contigent on args:

    # extract annotation files. Handle cases where either a parent dir or file is provided as an argument
    if args.no_annotation is not True:
        anno_file_subset = annotation_extract(args, gene_name, anno_file, gene_ids)

    # plot gene neighbourhood figure
    gff_input_file = gff_outpath
    if args.no_plot is not True:
        plot_gene_neighbourhood(args, gene_name, gff_input_file, goi_start, goi_end)

    # extract fasta files. Handle cases where either a parent dir or file is provided as an argument
    if args.no_fasta is not True:
        # NOTE: adding the gene IDs list without maturation genes
        fasta_file_subset = fasta_neighbourhood_extract(
            args,
            gene_name,
            fasta_file,
            anno_file_subset,
            gene_ids,
            gene_ids_nomaturation,
            gene_ids_ppi_candidates,
        )

    # call nife ssu extractor
    if args.no_ssu is not True:
        extract_nife_ssu(args, gene_name, gff_subset, fasta_file_subset)


def extract_nife_ssu(
    args: argparse.Namespace,
    target_name: str,
    gff_input: pl.DataFrame,
    fasta_file: Path,
) -> None:
    """Find the nearest NiFe SSU based on gene annotations in gff dataframe"""
    # params
    output_dir = args.output_dir
    gff = gff_input
    gene_name = target_name
    lsu_id = target_name
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

    goi = gff.filter(pl.col("attributes").str.contains(f"ID={lsu_id}"))

    scaffold = goi.item(0, "seqid")
    strand = goi.item(0, "strand")
    start = goi.item(0, "start")
    end = goi.item(0, "end")

    window_start = max(start - upstream_window, 0)
    window_end = end + downstream_window

    # extract neighbourhood
    neighbours = gff.filter(
        (pl.col("seqid") == scaffold)
        & (pl.col("start") <= window_end)
        & (pl.col("end") >= window_start)
        & (pl.col("strand") == strand)
    )

    # identify candidate matching codes
    ssu_candidates = neighbours.filter(
        pl.any_horizontal(
            [pl.col("attributes").str.contains(code) for code in NIFE_SSU_CODES]
        )
    )

    if not ssu_candidates.is_empty():
        # choose the closest nife ssu neighbour (i.e. smallest distance to target)
        ssu_candidates = ssu_candidates.with_columns(
            pl.when(pl.col("start") > end)
            .then(pl.col("start") - end)
            .otherwise(start - pl.col("end"))
            .alias("distance")
        )
        closest_ssu = ssu_candidates.sort("distance").head(1)

        # extract gene ID from attributes col
        ssu_id = (
            closest_ssu.select(
                pl.col("attributes").str.extract_groups(r"ID=([^;]+)").struct.field("1")
            )
            .to_series()
            .to_list()[0]
        )
        # make output path
        outpath1 = Path(output_dir) / gene_name / f"{ssu_id}-NiFe_SSU.faa"
        outpath1.parent.mkdir(exist_ok=True, parents=True)
        # make other output path
        outpath2 = Path(output_dir) / gene_name / f"{lsu_id}-{ssu_id}-NiFe_LSU_SSU.faa"
        outpath2.parent.mkdir(exist_ok=True, parents=True)

        # find subset of fasta file based on gene neighbourhood ids
        with open_fasta(fasta_file) as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            ssu_record = next((rec for rec in records if rec.id == ssu_id), None)
            lsu_record = next((rec for rec in records if rec.id == lsu_id), None)
            lsu_and_ssu_records = [lsu_record, ssu_record]
            # all_fastas = SeqIO.parse(handle, "fasta")
            # gene_neighbours = [record for record in all_fastas if record.id in gene_ids]
            with open(outpath1, "w") as f:
                SeqIO.write(ssu_record, f, "fasta")
            with open(outpath2, "w") as f:
                SeqIO.write(lsu_and_ssu_records, f, "fasta")

        # TODO: add boltz, chai, colabfold fastas of LSU+SSU
    else:
        print(f"No neighbours near {lsu_id} match any NiFe SSU annotation codes")


def plot_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_input_file: Path,
    goi_start: str,
    goi_end: str,
) -> None:
    """Plot the genetic neighbourhood around gene of interest using dna_features_viewer"""
    # params
    output_dir = args.output_dir
    dpi = args.plot_dpi
    format = args.plot_format
    target_gene_id = gene_name.split("___")[1]
    gff_input = str(gff_input_file)
    # gff_input = gff_input_file
    window_start = args.upstream_window
    window_end = args.downstream_window

    # Define custom class for dna_features_viewer. see: https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/index.html#custom-biopython-translators
    # TODO: update this with more colours/options
    class CustomTranslator(BiopythonTranslator):
        def compute_feature_color(self, feature):
            # Color the target gene green, annotation matches blue, and all others grey
            if "ID" in feature.qualifiers:
                gene_id = feature.qualifiers["ID"][0]
                gene_id = gene_id.split("___")[1]
                if "Name" in feature.qualifiers:
                    gene_anno = feature.qualifiers["Name"][0]
                    # colour genes that match annotations of familiar neighbours blue
                    if gene_anno in NIFE_GROUP1_PPI_CODES:
                        return "#cba6f7"
                    if gene_anno in NIFE_FRHB_CODES:
                        return "#b4befe"
                    if gene_anno in NIFE_SSU_CODES:
                        return "#89b4fa"
                    if gene_anno in NIFE_MATURATION_CODES:
                        return "#f9e2af"
                    # colour expected target with correct annotation as green
                    if gene_anno in NIFE_LSU_CODES:
                        return "#a6e3a1"
                if "product" in feature.qualifiers:
                    gene_desc = feature.qualifiers["product"][0]
                    if "maturation" in gene_desc:
                        return "#f9e2af"
                    if "chaperone" in gene_desc:
                        return "#f9e2af"
                # colour target gene red if it doesn't match anything else
                if target_gene_id == gene_id:
                    return "#f38ba8"
            return "#cdd6f4"

    # plot genetic neighbourhood figure, using CustomTranslator
    record = CustomTranslator().translate_record(
        record=gff_input, record_class="linear", filetype="gff"
    )

    # Suppose goi_start/goi_end are gene positions on this scaffold
    window_start = max(0, goi_start - args.upstream_window)
    window_end = goi_end + args.downstream_window
    # Clamp to sequence length to avoid overshooting
    window_end = min(window_end, record.sequence_length)

    # set sequence length of record to avoid out of bounds error
    record.sequence_length = max(window_end, record.sequence_length)
    # crop plot to window
    record = record.crop((window_start, window_end))

    # plot figure
    fig, ax = plt.subplots()
    ax, _ = record.plot(figure_width=20, strand_in_label_threshold=3)
    ax.figure.tight_layout()

    # save figure
    outpath = Path(output_dir) / gene_name / f"{gene_name}___plot.{format}"
    outpath.parent.mkdir(exist_ok=True, parents=True)
    ax.figure.savefig(outpath, dpi=dpi, format=format)
    # close figure to avoid memory issues
    plt.close("all")


def fasta_neighbourhood_extract(
    args: argparse.Namespace,
    gene_name: str,
    fasta_file: Path,
    annotation_file: Path,
    gene_ids_full: list,
    gene_ids_nomaturation: list,
    gene_ids_ppi_candidates: list,
) -> Path:
    """Extract fasta seqs surrounding gene of interest and output .faa files in a new dir, each containing a pairwise combination of fastas for AlphaPulldown input"""
    # params
    output_dir = args.output_dir
    fasta_file = fasta_file
    annotation_file = annotation_file
    gene_ids = gene_ids_full
    gene_ids_nomaturation = gene_ids_nomaturation
    gene_ids_ppi_candidates = gene_ids_ppi_candidates
    gene_name = gene_name

    # map annotation descriptions to gene ids
    df = pl.read_csv(annotation_file, separator="\t", has_header=True, quote_char=None)
    df = df.select(pl.col("gene_callers_id"), pl.col("function"))
    df = df.unique(subset=["gene_callers_id"], keep="first")
    desc_map = df.transpose(column_names="gene_callers_id").to_dict(as_series=False)

    # output all fastas from the entire gene neighbourhood
    outpath_all = Path(output_dir) / gene_name / f"{gene_name}-gene_neighbours_all.faa"
    outpath_all.parent.mkdir(exist_ok=True, parents=True)
    with open_fasta(fasta_file) as handle:
        records = list(SeqIO.parse(handle, "fasta"))
        fasta_subset = [rec for rec in records if rec.id in gene_ids]
        descriptive_fastas = [
            SeqRecord(
                seq=rec.seq,
                name=rec.id,
                id=rec.id,
                description=str(desc_map.get(rec.id, [None])[0]),
            )
            for rec in fasta_subset
        ]
        with open(outpath_all, "w") as f:
            SeqIO.write(descriptive_fastas, f, "fasta")

    # output all fastas minus the maturation factors (set1), use fasta file generated above as input
    if args.pairwise_set1 is True:
        outpath_set1 = (
            Path(output_dir) / gene_name / f"{gene_name}-gene_neighbours_set1.faa"
        )
        outpath_set1.parent.mkdir(exist_ok=True, parents=True)
        with open_fasta(outpath_all) as handle:
            fastafile = SeqIO.parse(handle, "fasta")
            fasta_subset = [rec for rec in fastafile if rec.id in gene_ids_nomaturation]
            with open(outpath_set1, "w") as f:
                SeqIO.write(fasta_subset, f, "fasta")
        # generate fasta pairwise combinations:
        pairwise_fasta_generation(args, gene_name, outpath_set1, "Set1")

    # output just the fastas that are most likley to be NiFe PPI candidates (set2)
    if args.pairwise_set2 is True:
        outpath_set2 = (
            Path(output_dir) / gene_name / f"{gene_name}-gene_neighbours_set2.faa"
        )
        outpath_set2.parent.mkdir(exist_ok=True, parents=True)
        with open_fasta(outpath_all) as handle:
            fastafile = SeqIO.parse(handle, "fasta")
            fasta_subset = [
                rec for rec in fastafile if rec.id in gene_ids_ppi_candidates
            ]
            with open(outpath_set2, "w") as f:
                SeqIO.write(fasta_subset, f, "fasta")
        # generate fasta pairwise combinations:
        pairwise_fasta_generation(args, gene_name, outpath_set2, "Set2")

    return outpath_all


def pairwise_fasta_generation(
    args: argparse.Namespace,
    target_id: str,
    input_fasta_set: Path,
    set_num: str,
) -> None:
    """Generate multiple fasta files containing pairwise combinations of all gene gene_neighbours"""
    # params
    output_dir = args.output_dir

    # Read all sequences into a list (records)
    records = list(SeqIO.parse(input_fasta_set, "fasta"))

    # generate fasta pairs using itertools
    pairs = combinations_with_replacement(records, 2)
    for i, (rec1, rec2) in enumerate(pairs, start=1):
        # Check if target is in this pair
        assert isinstance(rec1.name, str)
        assert isinstance(rec2.name, str)

        # get gene names to write new file names
        genome = rec1.name.split("___")[0]
        gene1 = rec1.name.split("___")[1]
        gene2 = rec2.name.split("___")[1]

        # Build output file name
        faaname = f"{genome}___{gene1}-{gene2}.faa"
        stemname = f"{genome}___{gene1}-{gene2}"

        # write pairwise fastas in chai format
        if args.chai_fastas is True:
            chai_records = [
                SeqRecord(
                    seq=record.seq,
                    name=record.id,
                    id=f"protein|{record.id}",
                    description="",
                )
                for record in records
            ]
            outpath = (
                Path(output_dir)
                / target_id
                / f"Chai1_fasta_pairs-{set_num}-{target_id}"
                / faaname
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write(chai_records, outpath, "fasta")

        # write pairwise fastas in boltz format
        if args.boltz_fastas is True:
            boltz_records = []
            for chain_id, rec in zip(ascii_uppercase, (rec1, rec2)):
                boltz_records.append(
                    SeqRecord(
                        seq=rec.seq,
                        name=rec.id,
                        id=f"{chain_id}|protein|MSAPATH/{stemname}.a3m",
                        description="",
                    )
                )
            outpath = (
                Path(output_dir)
                / target_id
                / f"Boltz2_fasta_pairs-{set_num}-{target_id}"
                / faaname
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write(boltz_records, outpath, "fasta")

        # write pairwise fastas in colabfold format
        if args.colabfold_fastas is True:
            outpath = (
                Path(output_dir)
                / target_id
                / f"ColabFold_fasta_pairs-{set_num}-{target_id}"
                / faaname
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            with open(outpath, "w") as f:
                f.write(f">{stemname}\n")
                f.write(f"{str(rec1.seq)}:\n")
                f.write(f"{str(rec2.seq)}")

        # write pairwise fastas in a regular format if no special formatting options are set, or add these too if --regular_fastas is set
        if (
            args.boltz_fastas and args.chai_fastas and args.colabfold_fastas is not True
        ) or args.regular_fastas is True:
            outpath = (
                Path(output_dir)
                / target_id
                / f"Fasta_pairs-{set_num}-{target_id}"
                / faaname
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write([rec1, rec2], outpath, "fasta")


def annotation_extract(
    args: argparse.Namespace,
    gene_name: str,
    annotation_file: Path,
    gene_ids: list,
) -> Path:
    """Extract annotation info surrounding gene of interest and output to new .tsv"""
    # params
    output_dir = args.output_dir
    gene_ids = gene_ids
    anno_file = annotation_file
    gene_name = gene_name

    # scan annotation tsv file
    anno = pl.scan_csv(anno_file, separator="\t", has_header=True, quote_char=None)

    # create subset of annotation table based on gene_ids list, and sorted based on gene_ids and evalues
    anno_subset = anno.filter(pl.col("gene_callers_id").is_in(gene_ids)).sort(
        by=["gene_callers_id", "e_value"], descending=[False, False]
    )
    # anno_subset = anno.filter(pl.col("gene_callers_id").is_in(gene_ids))

    # Save annotation tsv subset
    outpath = Path(output_dir) / gene_name / f"{gene_name}___annotations.tsv"
    outpath.parent.mkdir(exist_ok=True, parents=True)
    anno_subset.collect().write_csv(outpath, separator="\t", include_header=True)

    return outpath


def open_fasta(file_path: Path) -> TextIOWrapper:
    """Open a fasta file that might be gzipped or plain text"""
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")


def build_file_index(root_dir: Path) -> dict[str, list[Path]]:
    """Build index of mapping filenames to filepaths"""
    index = defaultdict(list)
    for path in root_dir.rglob("*"):
        if path.is_file():
            index[path.name].append(path.resolve())
    return index


def find_gff_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """Use prebuilt index to find gff file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}_cog.gff"
        or filename == f"{target_name}_cog.gff.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No file found for {target_name}.")
    print(f"gff file found: {matches[0]}")
    return matches[0]


def find_annotation_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """Use prebuilt index to find annotation file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}_annotations.tsv"
        or filename == f"{target_name}_annotations.tsv.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No file found for {target_name}.")
    print(f"annotation file found: {matches[0]}")
    return matches[0]


def find_fasta_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """Use prebuilt index to find fasta file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}.faa" or filename == f"{target_name}.faa.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No file found for {target_name}.")
    print(f"fasta file found: {matches[0]}")
    return matches[0]


def process_target_genes(args: argparse.Namespace) -> None:
    """Handles processing input depending on what arguments were parsed"""
    # Read or write index file of GlobDB
    if args.index_path:
        pickle_file = args.index_path
        if pickle_file.is_file():
            print("GlobDB index exists, reading...")
            with open(pickle_file, "rb") as f:
                file_index = pickle.load(f)
    else:
        print("Writing index from GlobDB dir provided. May take a while...")
        new_pickle_file = "./globdb_index.pkl"
        file_index = build_file_index(args.data_dir)
        with open(new_pickle_file, "wb") as f:
            pickle.dump(file_index, f)
            print(f"Index written to {new_pickle_file} in current directory")

    # Process target genes from target list
    if args.gene_list:
        with open(args.gene_list) as target_file:
            for line in target_file:
                gene_name = line.rstrip()
                print(f"Searching for {gene_name}")
                target_name = gene_name.split("___")[0]
                gff_file = find_gff_file_from_index(file_index, target_name)
                anno_file = find_annotation_file_from_index(file_index, target_name)
                fasta_file = find_fasta_file_from_index(file_index, target_name)
                extract_gene_neighbourhood(
                    args, gene_name, gff_file, anno_file, fasta_file
                )
    # Or process single target
    if args.gene_name:
        gene_name = args.gene_name.rstrip()
        print(f"Searching for {gene_name}")
        target_name = gene_name.split("___")[0]
        gff_file = find_gff_file_from_index(file_index, target_name)
        anno_file = find_annotation_file_from_index(file_index, target_name)
        fasta_file = find_fasta_file_from_index(file_index, target_name)
        extract_gene_neighbourhood(args, gene_name, gff_file, anno_file, fasta_file)


def main():
    """Handles broad control flow of all functions"""
    args = parse_arguments()
    process_target_genes(args)


if __name__ == "__main__":
    main()
