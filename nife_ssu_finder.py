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
from re import sub
from constants.nife_hmm_codes import (
    NIFE_LSU_CODES,
    NIFE_SSU_CODES,
    NIFE_FRHB_CODES,
    NIFE_GROUP2_PPI_CODES,
    NIFE_GROUP1_PPI_CODES,
    NIFE_GROUP3_PPI_CODES,
    NIFE_GROUP4_PPI_CODES,
    NIFE_MATURATION_CODES,
    NIFE_PPI_CANDIDATES,
    IRON_HYDROGENASE_CODES,
)
import shutil
from io import TextIOWrapper
from string import ascii_uppercase
import pickle
import gzip
import logging
import argparse
import textwrap
import polars as pl
from BCBio import GFF
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.font_manager as fm
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
from itertools import combinations_with_replacement


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
        help="Path to list of genes of interest, one per line (Either -l or -n is required).",
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
        "--add_hyddb",
        dest="add_hyddb",
        type=Path,
        required=False,
        metavar="FILE.tsv",
        help="Path to .tsv of gene ID's and respective HydDB classification",
    )

    parser.add_argument(
        "--add_arial_font",
        dest="add_arial_font",
        type=Path,
        required=False,
        metavar="FONT.ttf",
        help="Path to Arial font .ttf file to incorporate into final plot [Default: none]",
    )

    parser.add_argument(
        "-U",
        "--upstream",
        dest="upstream_window",
        type=int,
        default=10000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) upstream of GOI [Default: 5000]",
    )

    parser.add_argument(
        "-D",
        "--downstream",
        dest="downstream_window",
        type=int,
        default=10000,
        required=False,
        metavar="INT",
        help="Size of window (in base pairs) downstream of GOI [Default: 5000]",
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
        "-S",
        "--one_strand",
        dest="one_strand",
        action="store_true",
        help="Option to only extract genetic neighbourhood from the same strand the gene-of-interest is located on. Default is to extract neighbourhood of both strands [Default: off]",
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
        "--retain_full_gff",
        dest="retain_full_gff",
        action="store_true",
        help="Keep all gene annotations in .gff file, will plot all info too [Default: retain only best scoring evalue annotations in output .gff file]",
    )

    anno_args.add_argument(
        "--kofam_only",
        dest="kofam_only",
        action="store_true",
        help="Use only KOfam annotations in final gene neighbourhood plot [Default: uses best out of all types of annotations per gene]",
    )

    anno_args.add_argument(
        "--pfam_only",
        dest="pfam_only",
        action="store_true",
        help="Use only Pfam annotations in final gene neighbourhood plot [Default: uses best out of all types of annotations per gene]",
    )

    anno_args.add_argument(
        "--cogfun_only",
        dest="cogfun_only",
        action="store_true",
        help="Use only COG20_FUNCTION annotations in final gene neighbourhood plot [Default: uses best out of all types of annotations per gene]",
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

    return args


def extract_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_file: Path,
    anno_file: Path,
    fasta_file: Path,
    taxonomy_df: pl.DataFrame,
) -> None:
    """Finds the genes upstream and downstream of gene of interest in .gff file and returns/outputs new subset of the .gff file"""
    # params
    output_dir = args.output_dir  # Your output file
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

    # step 1: get gff subset
    raw_gff = pl.read_csv(
        gff_file,
        separator="\t",
        quote_char=None,
        comment_prefix="#",
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
    goi_row = raw_gff.filter(
        pl.col("attributes").str.contains(f"ID={gene_name}(?:;|$)")
    )

    # Get genetic neighbourhood coordinates to center on gene of interest
    goi_start = goi_row[0, "start"]
    goi_end = goi_row[0, "end"]
    goi_scaffold = goi_row[0, "seqid"]
    goi_strand = goi_row[0, "strand"]

    # Define window coordinates
    window_start: int = max(goi_start - upstream_window, 0)
    window_end: int = goi_end + downstream_window

    # Get subset of gff file based on window. Get genes from both strands, unless -S|--one_strand flag is provided
    if args.one_strand is True:
        # gff of all genes in the target neighbourhood
        gff_subset = raw_gff.filter(
            (pl.col("strand") == goi_strand)
            & (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
            & (pl.col("type") == "CDS")
        )
    else:
        # gff of all genes in the target neighbourhood
        gff_subset = raw_gff.filter(
            (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
            & (pl.col("type") == "CDS")
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
        .unique()
        .to_list()
    )

    # remove other attributes info and add key column to merge and sort on
    blank_gff = gff_subset.with_columns(
        pl.col("attributes").str.extract(r"ID=([^;]+)(?:;|$)", 1).alias("ID"),
        pl.col("attributes")
        .str.extract(r"ID=.*___(\d+)(?:;|$)", 1)
        .cast(pl.Int32)
        .alias("ID_num"),
    )

    # read in dataframe of annotation file subset
    anno_file_subset = annotation_extract(args, gene_name, anno_file, gene_ids)
    anno_df = pl.read_csv(
        anno_file_subset,
        separator="\t",
        quote_char=None,
        has_header=True,
        new_columns=["ID", "db_xref", "Name", "product", "evalue"],
        schema_overrides={"evalue": pl.Float64},
    )

    # clean rows that contain "!!!" in their names and description, only keeping the first name and description
    anno_df = anno_df.with_columns(
        pl.col("Name").str.extract(r"([^!]+)"),
        pl.col("product").str.extract(r"([^!]+)"),
    )

    # filter to use only certain HMM codes
    anno_df = anno_df  # default
    if args.pfam_only:
        anno_df = anno_df.filter(pl.col("db_xref") == "Pfam")
    if args.kofam_only:
        anno_df = anno_df.filter(pl.col("db_xref") == "KOfam")
    if args.cogfun_only:
        anno_df = anno_df.filter(pl.col("db_xref") == "COG20_FUNCTION")

    # merge gff_subset with anno_df on shared gene ID
    merge_gff_anno = blank_gff.join(anno_df, on="ID", how="left")

    rebuilt_gff_full = (
        merge_gff_anno.with_columns(
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
        .select(
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
                "evalue",
                "ID",
                "ID_num",
            ]
        )
        .sort(by=["ID_num", "evalue"])
    )

    # keep only the best annotation row per ID
    rebuilt_gff_uniq = (
        rebuilt_gff_full.unique(subset="ID", keep="first")
        .unique(subset="ID", keep="first")
        .sort(by="ID_num")
        .select(
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
    )

    # TODO: figure out what to do with this...
    # get a gene list of just SSU candidates
    gff_subset_ssus = rebuilt_gff_full.filter(
        pl.col("attributes").str.contains_any(NIFE_SSU_CODES)
        & (~pl.col("attributes").str.contains(rf"ID={gene_name}(?:;|$)"))
    )
    # gene_ids_ssu_candidates: list = (
    #     gff_subset_ssus.select(
    #         pl.col("attributes")
    #         .str.extract_groups(r"ID=([^;]+)")
    #         .struct.field("1")
    #         .alias("gene_id")
    #     )
    #     .to_series()
    #     .unique()
    #     .to_list()
    # )

    # gff subset WITHOUT maturation genes, for smaller fasta file extraction
    gff_subset_nomaturation = rebuilt_gff_full.filter(
        (~pl.col("attributes").str.contains_any(NIFE_MATURATION_CODES))
    )
    gene_ids_nomaturation: list = (
        gff_subset_nomaturation.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .unique()
        .to_list()
    )

    # gff subset WITH ONLY relevant PPI candidates for even smaller fasta file candidates
    gff_subset_ppi_candidates = rebuilt_gff_full.filter(
        (pl.col("attributes").str.contains_any(NIFE_PPI_CANDIDATES))
    )
    gene_ids_ppi_candidates: list = (
        gff_subset_ppi_candidates.select(
            pl.col("attributes")
            .str.extract_groups(r"ID=([^;]+)")
            .struct.field("1")
            .alias("gene_id")
        )
        .to_series()
        .unique()
        .to_list()
    )

    # Write output
    gff_dataframe_to_write = rebuilt_gff_uniq
    if args.kofam_only:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}-geneneighbours_with_KOfam_annotations.gff"
        )
    elif args.pfam_only:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}-geneneighbours_with_Pfam_annotations.gff"
        )
    elif args.cogfun_only:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}-geneneighbours_with_COGFUN_annotations.gff"
        )
    elif args.retain_full_gff:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}-geneneighbours_with_all_annotations.gff"
        )
        gff_dataframe_to_write = (
            rebuilt_gff_full  # keep all non-unique rows in gff file
        )
    else:
        gff_outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}-geneneighbours_with_best_annotations.gff"
        )

    gff_outpath.parent.mkdir(exist_ok=True, parents=True)

    gff_dataframe_to_write.write_csv(
        gff_outpath,
        separator="\t",
        include_header=False,
    )
    # TODO: make a logging file
    print(
        f"Extracted region ({window_start}-{window_end}) on {goi_scaffold} written to: {gff_outpath.name}"
    )

    # Execute downstream functions:

    # plot gene neighbourhood figure
    gff_input_file = gff_outpath
    if args.no_plot is not True:
        plot_gene_neighbourhood(
            args, gene_name, gff_input_file, goi_start, goi_end, taxonomy_df
        )

    # extract fasta sequences
    fasta_file_subset = fasta_neighbourhood_extract(
        args,
        gene_name,
        fasta_file,
        anno_file_subset,
        gene_ids,
        gene_ids_nomaturation,
        gene_ids_ppi_candidates,
    )

    # call nife ssu extractor, uses fasta file generated from previous step as input
    if args.no_ssu is not True:
        extract_nife_ssu(args, gene_name, rebuilt_gff_full, fasta_file_subset)


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

    goi = gff.filter(pl.col("attributes").str.contains(rf"ID={lsu_id}(?:;|$)"))

    scaffold = goi.item(0, "seqid")
    start = goi.item(0, "start")
    end = goi.item(0, "end")
    # strand = goi.item(0, "strand")

    window_start = max(start - upstream_window, 0)
    window_end = end + downstream_window

    # extract neighbourhood
    neighbours = gff.filter(
        (pl.col("seqid") == scaffold)
        & (pl.col("start") <= window_end)
        & (pl.col("end") >= window_start)
        & (pl.col("type") == "CDS")
        & (~pl.col("attributes").str.contains(rf"ID={lsu_id}(?:;|$)"))
        # & (pl.col("strand") == strand)
    )

    # identify candidate matching codes
    # TODO: use annotation table instead of single attributes column... do same for gff_extract function
    ssu_candidates = neighbours.filter(
        pl.any_horizontal(
            [pl.col("attributes").str.contains(code) for code in NIFE_SSU_CODES]
        )
    )

    # make summary table of LSU + SSU info
    summary_table = Path(output_dir) / "NiFe_LSU_SSU_pairs.tsv"

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

        # write to summary_table
        with open(summary_table, "a") as f:
            f.write(f"{lsu_id}\t{ssu_id}\n")

        # set output file names
        genome_name = lsu_id.split("___")[0]
        lsu_id_num = lsu_id.split("___")[1]
        ssu_id_num = ssu_id.split("___")[1]
        lsu_ssu_faaname = f"{genome_name}___{lsu_id_num}-{ssu_id_num}"
        lsu_ssu_faaname2 = f"{genome_name}___{ssu_id_num}-{lsu_id_num}"
        lsu_ssu_file_names = [lsu_ssu_faaname, lsu_ssu_faaname2]

        # save NiFe SSU fasta file
        outpath = Path(output_dir) / gene_name / f"{ssu_id}-NiFe_SSU.faa"
        outpath.parent.mkdir(exist_ok=True, parents=True)
        with open_fasta(fasta_file) as handle:
            records = list(SeqIO.parse(handle, "fasta"))
            ssu_record = next((rec for rec in records if rec.id == ssu_id), None)
            with open(outpath, "w") as f:
                SeqIO.write(ssu_record, f, "fasta")

        # if pairwise_set1 or pairwise_set2 options are set, then LSU-SSU fastas already exist and just needs to be copied to a new dir
        if args.pairwise_set1 or args.pairwise_set2 is True:
            outpath = Path(output_dir) / gene_name / "NiFe_LSU_SSU_pair"
            outpath.mkdir(exist_ok=True, parents=True)
            target_path = Path(output_dir) / gene_name
            for file in target_path.rglob("*"):
                if outpath in file.parents:
                    continue
                file_name = file.stem.split("_AFformat_")[0]
                if file.is_file() and file_name in lsu_ssu_file_names:
                    shutil.copy2(file, outpath)
        else:
            # write NiFe_LSU_SSU fasta file from scratch if pairwise_set options are not set
            outpath2 = (
                Path(output_dir) / gene_name / f"{lsu_ssu_faaname}-NiFe_LSU_SSU.faa"
            )
            outpath2.parent.mkdir(exist_ok=True, parents=True)
            with open_fasta(fasta_file) as handle:
                records = list(SeqIO.parse(handle, "fasta"))
                ssu_record = next((rec for rec in records if rec.id == ssu_id), None)
                lsu_record = next((rec for rec in records if rec.id == lsu_id), None)
                lsu_and_ssu_records = [lsu_record, ssu_record]
                with open(outpath2, "w") as f:
                    SeqIO.write(lsu_and_ssu_records, f, "fasta")

    else:
        # TODO: make a logging file
        print(f"No neighbours near {lsu_id} match any NiFe SSU annotation codes")

        # write to summary_table
        with open(summary_table, "a") as f:
            f.write(f"{lsu_id}\t-\n")


def plot_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_input_file: Path,
    goi_start: str,
    goi_end: str,
    taxonomy_df: pl.DataFrame,
) -> None:
    """Plot the genetic neighbourhood around gene of interest using dna_features_viewer"""
    # params
    output_dir = args.output_dir
    dpi = args.plot_dpi
    format = args.plot_format
    target_gene_id = gene_name.split("___")[1]
    genome_name = gene_name.split("___")[0]
    gff_input = str(gff_input_file)
    window_start = args.upstream_window
    window_end = args.downstream_window

    # add font if provided:
    if args.add_arial_font:
        arial_font = args.add_arial_font
        # arial_font_bold = "/home/james/Downloads/Arial Bold.ttf"
        fm.fontManager.addfont(arial_font)
        # fm.fontManager.addfont(arial_font_bold)
        rcParams["font.sans-serif"] = "Arial"
        rcParams["font.family"] = "Arial"
        rcParams["font.size"] = 10

    # Define custom class for dna_features_viewer. see: https://edinburgh-genome-foundry.github.io/DnaFeaturesViewer/index.html#custom-biopython-translators
    # TODO: update this with more colours/options
    class CustomTranslator(BiopythonTranslator):
        def compute_feature_color(self, feature):
            if "ID" in feature.qualifiers:
                gene_id = feature.qualifiers["ID"][0].split("___")[1]
                if "product" in feature.qualifiers:
                    gene_desc = feature.qualifiers["product"][0]
                    if "maturation" in gene_desc:
                        return "#f9e2af"
                    if "chaperone" in gene_desc:
                        return "#f9e2af"
                if "Name" in feature.qualifiers:
                    gene_anno = feature.qualifiers["Name"][0]
                    # colour genes that match annotations of familiar neighbours
                    if gene_anno in IRON_HYDROGENASE_CODES:
                        return "#f5a97f"
                    if gene_anno in NIFE_GROUP1_PPI_CODES:
                        return "#cba6f7"
                    if gene_anno in NIFE_GROUP2_PPI_CODES:
                        return "#f5c2e7"
                    if gene_anno in NIFE_GROUP3_PPI_CODES:
                        return "#94e2d5"
                    if gene_anno in NIFE_GROUP4_PPI_CODES:
                        return "#7287fd"
                    if gene_anno in NIFE_FRHB_CODES:
                        return "#94e2d5"
                    if gene_anno in NIFE_SSU_CODES:
                        return "#89b4fa"
                    if gene_anno in NIFE_MATURATION_CODES:
                        return "#f9e2af"
                    # colour expected target with correct annotation as green
                    if gene_anno in NIFE_LSU_CODES:
                        return "#a6e3a1"
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

    # get taxonomy/species name of genome
    species = (
        taxonomy_df.filter(pl.col("genome_id") == f"{genome_name}")
        .select("species")
        .item()
    )

    figure_text = f"Species: {species}; Genome: {genome_name}; Gene: {target_gene_id}"

    # get hydrogenase classification
    if args.add_hyddb:
        hyddb_df = pl.read_csv(
            args.add_hyddb,
            separator="\t",
            has_header=False,
            new_columns=["gene_id", "classification"],
        )
        hydclass = (
            hyddb_df.filter(pl.col("gene_id") == f"{gene_name}")
            .get_column("classification")
            .first()
        )
        if hydclass is not None:
            figure_text = f"Species: {species}; Genome: {genome_name}; Gene: {target_gene_id}; [NiFe] group {hydclass}"

    # plot figure
    fig, ax = plt.subplots()
    ax, _ = record.plot(figure_width=20, strand_in_label_threshold=3)
    ax.figure.tight_layout()

    # add genome name to figure
    plt.figtext(
        0.02,
        -0.05,
        figure_text,
        ha="left",
        va="top",
        fontsize=12,
        style="normal",
    )

    # save figure
    outpath = Path(output_dir) / gene_name / f"{gene_name}-plot.{format}"
    outpath.parent.mkdir(exist_ok=True, parents=True)
    ax.figure.savefig(outpath, dpi=dpi, format=format, bbox_inches="tight")
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
                    seq=rec.seq,
                    name=rec.id,
                    id=f"protein|{rec.id}",
                    description="",
                )
                for rec in (rec1, rec2)
            ]
            outpath = (
                Path(output_dir)
                / target_id
                / f"Chai1_fasta_pairs-{set_num}-{target_id}"
                / f"{stemname}_AFformat_Chai1.faa"
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
                / f"{stemname}_AFformat_Boltz2.faa"
            )
            outpath.parent.mkdir(exist_ok=True, parents=True)
            SeqIO.write(boltz_records, outpath, "fasta")

        # write pairwise fastas in colabfold format
        if args.colabfold_fastas is True:
            outpath = (
                Path(output_dir)
                / target_id
                / f"ColabFold_fasta_pairs-{set_num}-{target_id}"
                / f"{stemname}_AFformat_ColabFold.faa"
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
    outpath = Path(output_dir) / gene_name / f"{gene_name}-annotations.tsv"
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
        raise FileNotFoundError(f"No gff file found for {target_name}.")
    return matches[0]


def find_annotation_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """use prebuilt index to find annotation file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}_annotations.tsv"
        or filename == f"{target_name}_annotations.tsv.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"No annotations file found for {target_name}.")
    return matches[0]


def find_fasta_file_from_index(
    file_index: dict[str, list[Path]], target_name: str
) -> Path:
    """use prebuilt index to find fasta file of interest"""
    matches = [
        path
        for filename, paths in file_index.items()
        if filename == f"{target_name}.faa" or filename == f"{target_name}.faa.gz"
        for path in paths
    ]
    if not matches:
        raise FileNotFoundError(f"no file found for {target_name}.")
    return matches[0]


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
        raise FileNotFoundError("ERROR: No taxonomy file found")
    return matches[0]


def process_target_genes(args: argparse.Namespace) -> None:
    """handles processing input depending on what arguments were parsed"""
    # read or write index file of globdb
    file_index = None
    if args.index_path:
        pickle_file = args.index_path
        if pickle_file.is_file():
            print("GlobDB index exists, reading...")
            with open(pickle_file, "rb") as f:
                file_index = pickle.load(f)
    else:
        print("writing index from globdb dir provided. may take a while...")
        new_pickle_file = "./globdb_index.pkl"
        file_index = build_file_index(args.data_dir)
        with open(new_pickle_file, "wb") as f:
            pickle.dump(file_index, f)
            print(f"index written to {new_pickle_file} in current directory")

    # TODO:
    # find and process taxonomy.tsv file once
    taxonomy_file = find_taxonomy_file_from_index(file_index)
    taxonomy_df = taxonomy_dataframe(taxonomy_file)

    # process target genes from target list
    if args.gene_list:
        with open(args.gene_list) as target_file:
            for line in target_file:
                gene_name = line.rstrip()
                print(f"searching {gene_name}")
                target_name = gene_name.split("___")[0]
                try:
                    gff_file = find_gff_file_from_index(file_index, target_name)
                except FileNotFoundError as e:
                    print(f"warning: {e}. skipping.")
                    continue
                anno_file = find_annotation_file_from_index(file_index, target_name)
                fasta_file = find_fasta_file_from_index(file_index, target_name)
                extract_gene_neighbourhood(
                    args, gene_name, gff_file, anno_file, fasta_file, taxonomy_df
                )

    # or process single target
    if args.gene_name:
        gene_name = args.gene_name.rstrip()
        print(f"searching {gene_name}")
        target_name = gene_name.split("___")[0]
        try:
            gff_file = find_gff_file_from_index(file_index, target_name)
        except FileNotFoundError as e:
            print(f"warning: {e} skipping.")
            return
        anno_file = find_annotation_file_from_index(file_index, target_name)
        fasta_file = find_fasta_file_from_index(file_index, target_name)
        extract_gene_neighbourhood(
            args, gene_name, gff_file, anno_file, fasta_file, taxonomy_df
        )


def main():
    """handles broad control flow of all functions"""
    args = parse_arguments()
    process_target_genes(args)


if __name__ == "__main__":
    main()
