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
import gzip
import logging
import argparse
import textwrap
import pandas as pd
import polars as pl
import seaborn as sns
from BCBio import GFF
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator
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
    input_args = parser.add_mutually_exclusive_group(required=True)
    anno_args = parser.add_mutually_exclusive_group()

    input_args.add_argument(
        "-l",
        "--gene-list",
        dest="gene_list",
        type=Path,
        metavar="FILE",
        help="Path to list of genes of interest (Either -l or -n are required).",
    )

    input_args.add_argument(
        "-n",
        "--gene-name",
        dest="gene_name",
        type=str,
        metavar="STR",
        help="Name of gene of interest (Either -l or -n are required)",
    )

    # Add individual arguments
    parser.add_argument(
        "-d",
        "--datadir",
        dest="data_dir",
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
        "--add-annotation",
        dest="add_annotation",
        action="store_true",
        help="Option to return a table of annotations for the gene neighbourhood [Default: off]",
    )

    parser.add_argument(
        "-f",
        "--add-fasta",
        dest="add_fasta",
        action="store_true",
        help="Path to dir containing fasta files",
    )

    parser.add_argument(
        "-I",
        "--dpi",
        dest="plot_dpi",
        type=int,
        default=300,
        required=False,
        metavar="INT",
        help="Resolution of .png output plot in dpi [Default: 300]",
    )

    parser.add_argument(
        "-F",
        "--plot-format",
        dest="plot_format",
        choices=["png", "pdf", "svg"],
        default="png",
        required=False,
        metavar="STR",
        help="File format for output plot. Choices = png, pdf, svg [Default: png]",
    )

    parser.add_argument(
        "-N",
        "--no-plot",
        dest="no_plot",
        action="store_true",
        help="Option to not output genetic neighbourhood plots [Default: off]",
    )

    parser.add_argument(
        "-B",
        "--both-strands",
        dest="both_strands",
        action="store_true",
        help="Option to include both strands of genetic neighbourhood in output, rather than only genetic neighbours on the same strand as the gene-of-interest. [Default: off]",
    )

    parser.add_argument(
        "-p",
        "--add-pairwise-fastas",
        dest="add_pairwise_fastas",
        action="store_true",
        help="Option to generate fasta files containing pairwise combinations of all sequences in gene neighbourhood [Default: off]",
    )

    parser.add_argument(
        "-c",
        "--chai-fastas",
        dest="chai_fastas",
        action="store_true",
        help="Option to prepend fasta ID headers with '>protein|' for Chai-1 input [Default: off]",
    )

    anno_args.add_argument(
        "-K",
        "--use-kofam-annotation",
        dest="use_kofam_annotation",
        action="store_true",
        help="Replace default COG20_FUNCTION annotations with KOfam annotations",
    )

    anno_args.add_argument(
        "-P",
        "--use-pfam-annotation",
        dest="use_pfam_annotation",
        action="store_true",
        help="Replace default COG20_FUNCTION annotations with Pfam annotations",
    )

    # Parse arguments into args object
    args = parser.parse_args()

    # Info about input files
    if args.data_dir.is_file():
        print(f"Input .gff file provided: {args.data_dir}")

    if args.data_dir.is_dir():
        print(f"Target dir for .gff files provided. Searching {args.data_dir}/...")

    return args


def replace_gff_attributes(
    args: argparse.Namespace, gene_name: str, gff_file: Path, annotation_file: Path
) -> pl.DataFrame:
    """Rewrite the attributes column of the gff file with KOfam annotation information instead of the default COG20_FUNCTION annotations"""
    # params
    gff_file = gff_file
    annotation_file = annotation_file
    gene_name = gene_name
    output_dir = args.output_dir

    # load master annotation table of genome
    df = pl.read_csv(
        annotation_file,
        separator="\t",
        skip_lines=1,
        has_header=False,
        new_columns=["ID", "db_xref", "Name", "product", "evalue"],
        schema_overrides={"evalue": pl.Float64},
    )

    # keep only rows of KOfam or Pfam annotations, and only the best evalue per gene ID
    if args.use_kofam_annotation is True:
        df = df.filter(pl.col("db_xref") == "KOfam")
    if args.use_pfam_annotation is True:
        df = df.filter(pl.col("db_xref") == "Pfam")

    df = df.sort("evalue").group_by("ID").first()

    # load .gff file of genome
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
    args: argparse.Namespace, gene_name: str, gff_file: Path
) -> None:
    """Finds the genes upstream and downstream of gene of interest in .gff file and returns/outputs new subset of the .gff file"""
    # params
    gff_file = gff_file  # The input gff file
    output_dir = args.output_dir  # Your output file
    gene_name = gene_name
    genome_name = gene_name.split("___")[0]
    upstream_window = args.upstream_window  # Size of window upstream/downstream
    downstream_window = args.downstream_window  # Size of window upstream/downstream

    anno_file = find_annotation_file(args.data_dir, genome_name)

    # load gff dataframe
    if args.use_kofam_annotation is True:
        gff = replace_gff_attributes(args, gene_name, gff_file, anno_file)
    if args.use_pfam_annotation is True:
        gff = replace_gff_attributes(args, gene_name, gff_file, anno_file)
    else:
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

    # Get subset of gff file based on window. Get only genes from same strand as GOI, unless specified otherwise
    if args.both_strands is not True:
        gff_subset = gff.filter(
            (pl.col("strand") == goi_strand)
            & (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
        )
    else:
        gff_subset = gff.filter(
            (pl.col("seqid") == goi_scaffold)
            & (pl.col("start") <= window_end)
            & (pl.col("end") >= window_start)
        )
    # get subset of gff without maturation genes, for smaller fasta file extraction
    # TODO: incorporate this function in a smarter way:
    gff_subset_nomaturation = gff.filter(
        (pl.col("strand") == goi_strand)
        & (pl.col("seqid") == goi_scaffold)
        & (pl.col("start") <= window_end)
        & (pl.col("end") >= window_start)
        & (~pl.col("attributes").str.contains("maturation"))
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

    # Write output
    if args.use_kofam_annotation is True:
        outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_KOfam_annotations.gff"
        )
    if args.use_pfam_annotation is True:
        outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_Pfam_annotations.gff"
        )
    else:
        outpath = (
            Path(output_dir)
            / gene_name
            / f"{gene_name}___gene_neighbours_COG20_annotations.gff"
        )

    if outpath.exists():
        logging.warning(f"Writing over existing output .gff file: {outpath.name}")
    outpath.parent.mkdir(exist_ok=True, parents=True)

    gff_subset.write_csv(
        outpath,
        separator="\t",
        include_header=False,
    )
    print(
        f"Extracted region ({window_start}-{window_end}) on {goi_scaffold} written to: {outpath.name}"
    )

    # Execute downstream functions, contigent on args:

    # plot gene neighbourhood figure
    gff_input_file = str(outpath)
    if args.no_plot is not True:
        plot_gene_neighbourhood(
            args, gene_name, gff_input_file, window_start, window_end
        )

    # extract fasta files. Handle cases where either a parent dir or file is provided as an argument
    if args.add_fasta is True:
        target_name = gene_name.split("___")[0]
        fasta_file = find_fasta_file(args.data_dir, target_name)
        # NOTE: adding the gene IDs list without maturation genes
        fasta_neighbourhood_extract(
            args, gene_name, fasta_file, gene_ids, gene_ids_nomaturation
        )

    # extract annotation files. Handle cases where either a parent dir or file is provided as an argument
    if args.add_annotation is True:
        target_name = gene_name.split("___")[0]
        anno_file = find_annotation_file(args.data_dir, target_name)
        annotation_extract(args, gene_name, anno_file, gene_ids)


def plot_gene_neighbourhood(
    args: argparse.Namespace,
    gene_name: str,
    gff_input_file: str,
    window_start: int,
    window_end: int,
) -> None:
    """Plot the genetic neighbourhood around gene of interest using dna_features_viewer"""
    # params
    output_dir = args.output_dir
    dpi = args.plot_dpi
    format = args.plot_format
    target_gene_id = gene_name.split("___")[1]
    gff_input_file = gff_input_file
    window_start = window_start
    window_end = window_end
    # WARN: annotation code list is hardcoded... change that
    # TODO: provide a way to update this list from CLI args
    anno_codes_targets = [
        "COG0374",
        "COG3259",
        "COG4042",
        "COG3261",
    ]  # annotation codes associated with NiFe LSU
    anno_codes_neighbours_close = [
        "COG1035",
        "COG1740",
        "COG1908",
        "COG1941",
        "COG3260",
    ]  # annotation codes associated with NiFe SSU & FrhD/MvhD
    anno_codes_neighbours_other = [
        "COG5557",  # membrane subunit HybB
        "COG3658",  # cytochrome b
        "COG1969",  # HyaC cytochrome b
    ]

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
                    if gene_anno in anno_codes_neighbours_other:
                        return "#cba6f7"
                    if gene_anno in anno_codes_neighbours_close:
                        return "#89b4fa"
                    # colour expected target with correct annotation as green
                    if gene_anno in anno_codes_targets:
                        return "#a6e3a1"
                if "product" in feature.qualifiers:
                    gene_desc = feature.qualifiers["product"][0]
                    if "maturation" in gene_desc:
                        return "#f9e2af"
                # colour target gene red if it doesn't match anything else
                if target_gene_id == gene_id:
                    return "#f38ba8"
            return "#cdd6f4"

    # plot genetic neighbourhood figure, using CustomTranslator
    record = CustomTranslator().translate_record(gff_input_file)
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
    if outpath.exists():
        logging.warning(f"Writing over existing figure: {outpath.name}")
    outpath.parent.mkdir(exist_ok=True, parents=True)
    ax.figure.savefig(outpath, dpi=dpi, format=format)
    # close figure to avoid memory issues
    plt.close("all")


def fasta_neighbourhood_extract(
    args: argparse.Namespace,
    gene_name: str,
    fasta_file: Path,
    gene_ids_full: list,
    gene_ids_subset: list,
) -> None:
    """Extract fasta seqs surrounding gene of interest and output .faa files in a new dir, each containing a pairwise combination of fastas for AlphaPulldown input"""
    # params
    output_dir = args.output_dir
    fasta_file = fasta_file
    gene_ids = gene_ids_full
    gene_ids_subset = gene_ids_subset
    gene_name = gene_name

    # make output dir
    outpath = Path(output_dir) / gene_name / f"{gene_name}___gene_neighbours.faa"
    if outpath.exists():
        logging.warning(f"Writing over existing fastas: {outpath.name}")
    outpath.parent.mkdir(exist_ok=True, parents=True)

    # find subset of fasta file based on gene neighbourhood ids
    with open_fasta(fasta_file) as handle:
        all_fastas = SeqIO.parse(handle, "fasta")
        gene_neighbours = [record for record in all_fastas if record.id in gene_ids]
        with open(outpath, "w") as f:
            SeqIO.write(gene_neighbours, f, "fasta")

    # Do the same for subset of fasta list
    # make output dir
    outpath2 = (
        Path(output_dir) / gene_name / f"{gene_name}___gene_neighbours_subset.faa"
    )
    if outpath2.exists():
        logging.warning(f"Writing over existing fastas: {outpath2.name}")
    outpath2.parent.mkdir(exist_ok=True, parents=True)

    # find subset of fasta file based on gene neighbourhood ids
    with open_fasta(fasta_file) as handle:
        all_fastas = SeqIO.parse(handle, "fasta")
        gene_neighbours = [
            record for record in all_fastas if record.id in gene_ids_subset
        ]
        with open(outpath2, "w") as f:
            SeqIO.write(gene_neighbours, f, "fasta")

    # generate fasta pairwise combinations:
    # WARN: only doing it for subset of all gene neighbours
    # TODO: update this for better arg parsing
    if args.add_pairwise_fastas is True:
        fasta_pairwise_generation(args, outpath2, gene_name)


def fasta_pairwise_generation(
    args: argparse.Namespace, input_fasta: Path, goi_name: str
) -> None:
    """Generate multiple fasta files containing pairwise combinations of all gene gene_neighbours"""
    # params
    output_dir = args.output_dir
    target_id = goi_name
    input_fasta = input_fasta

    # Read all sequences into a list (records)
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # define local function for finding absolute difference between GOI and neighbouring genes based on GlobDB fasta ID number
    def dist_from_goi(record: SeqRecord):
        goi_num = int(target_id.split("___")[1])
        neighbour_num = int(record.id.split("___")[1])
        return abs(neighbour_num - goi_num)

    # sort all fasta records by the dist_from_goi function
    sorted_records = sorted(records, key=dist_from_goi)

    # Modify fasta header IDs of sorted_records for Chai-1 input, using the SeqRecord function
    if args.chai_fastas is True:
        seq_records = [
            SeqRecord(
                seq=record.seq,
                name=record.id,
                id=f"protein|{record.id}",
                description="",
            )
            for record in sorted_records
        ]
    else:
        seq_records = sorted_records

    # generate fasta pairs using itertools
    pairs = combinations_with_replacement(seq_records, 2)
    for i, (rec1, rec2) in enumerate(pairs, start=1):
        # Check if target is in this pair
        includes_target = target_id in {rec1.name, rec2.name}
        assert isinstance(rec1.name, str)
        assert isinstance(rec2.name, str)

        # get gene names to write new file names
        genome = rec1.name.split("___")[0]
        gene1 = rec1.name.split("___")[1]
        gene2 = rec2.name.split("___")[1]

        # Build output file name
        if includes_target:
            fname = f"{i:03d}_{genome}___{gene1}-{gene2}_geneofinterest.faa"
        else:
            fname = f"{i:03d}_{genome}___{gene1}-{gene2}.faa"

        # make output dir
        outpath = Path(output_dir) / target_id / f"Fasta_pairs_{target_id}" / fname
        if outpath.exists():
            logging.warning(f"Writing over existing fastas: {outpath.name}")
        outpath.parent.mkdir(exist_ok=True, parents=True)

        # Write the pairwise fasta file
        SeqIO.write([rec1, rec2], outpath, "fasta")


def annotation_extract(
    args: argparse.Namespace,
    gene_name: str,
    annotation_file: Path,
    gene_ids: list,
) -> None:
    """Extract annotation info surrounding gene of interest and output to new .tsv"""
    # params
    output_dir = args.output_dir
    gene_ids = gene_ids
    anno_file = annotation_file
    gene_name = gene_name
    anno = pd.read_csv(anno_file, delimiter="\t", compression="infer")

    # Reformat annotation table
    # create subset of annotation table based on gene_ids list, and sorted based on gene_ids and evalues
    anno_subset = anno[anno["gene_callers_id"].isin(gene_ids)].sort_values(
        by=["gene_callers_id", "e_value"], ascending=[True, True]
    )
    # make dataframe of gene_ids to merge with annotation table later
    gene_ids_df = pd.DataFrame({"gene_callers_id": gene_ids})
    # Merge gene_ids df with annotation subset df
    merged_with_missing = gene_ids_df.merge(
        anno_subset, on="gene_callers_id", how="left"
    )
    # Fill in missing rows from merge with NaN
    merged_filled = merged_with_missing.fillna(
        {
            "source": "-",
            "accession": "-",
            "function": "-",
            "e_value": "-",
        }
    )
    # Save annotation table
    outpath = Path(output_dir) / gene_name / f"{gene_name}___annotations.tsv"
    if outpath.exists():
        logging.warning(f"Writing over existing figure: {outpath.name}")
    outpath.parent.mkdir(exist_ok=True, parents=True)
    merged_filled.to_csv(outpath, sep="\t", index=False)


def find_gff_file(input_path: Path, target_name: str) -> Path:
    """Find target .gff file in input directory based on name of gene of interest"""
    matched_file = list(input_path.rglob(f"{target_name}*.gff*"))
    if not matched_file:
        raise FileNotFoundError(f"Could not find target file with name: {target_name}.")
    elif len(matched_file) > 1:
        logging.warning(f"Multiple files named: {target_name}: {matched_file}")
        logging.warning(f"Using first matched .gff file: {matched_file[0]}")
    return matched_file[0]


def find_fasta_file(input_path: Path, target_name: str) -> Path:
    """Find target .faa (or .fasta) file in input directory based on name of gene of interest"""
    matched_file = list(input_path.rglob(f"{target_name}*.faa*")) + list(
        input_path.rglob(f"{target_name}*.fasta*")
    )
    if not matched_file:
        raise FileNotFoundError(
            f"Could not find target .faa file with name: {target_name}."
        )
    elif len(matched_file) > 1:
        logging.warning(f"Multiple .faa files named: {target_name}: {matched_file}")
        logging.warning(f"Using first matched .faa file: {matched_file[0]}")
    return matched_file[0]


def open_fasta(file_path: Path) -> TextIOWrapper:
    """Open a fasta file that might be gzipped or plain text"""
    if file_path.suffix == ".gz":
        return gzip.open(file_path, "rt")
    else:
        return open(file_path, "r")


def find_annotation_file(input_path: Path, target_name: str) -> Path:
    """Find target annotation .tsv file in input directory based on name of gene of interest"""
    matched_file = list(input_path.rglob(f"{target_name}*.tsv*"))
    if not matched_file:
        raise FileNotFoundError(f"Could not find target file with name: {target_name}.")
    elif len(matched_file) > 1:
        logging.warning(
            f"Multiple annotation files named: {target_name}: {matched_file}"
        )
        logging.warning(f"Using first matched .tsv file: {matched_file[0]}")
    return matched_file[0]


def process_target_genes(args: argparse.Namespace) -> None:
    """Handles processing input depending on what arguments were parsed"""
    # Case 1: Passing list of genes of interest with general gff parent dir
    if args.gene_list and args.data_dir.is_dir():
        with open(args.gene_list) as target_file:
            for line in target_file:
                gene_name = line.rstrip()
                target_name = gene_name.split("___")[0]
                gff_file = find_gff_file(args.data_dir, target_name)
                extract_gene_neighbourhood(args, gene_name, gff_file)

    # Case 2: Passing singular gene of interest with general gff parent dir
    if args.gene_name and args.data_dir.is_dir():
        gene_name = args.gene_name
        target_name = gene_name.split("___")[0]
        gff_file = find_gff_file(args.data_dir, target_name)
        extract_gene_neighbourhood(args, gene_name, gff_file)


def main():
    """Handles broad control flow of all functions"""
    args = parse_arguments()
    process_target_genes(args)


if __name__ == "__main__":
    main()
