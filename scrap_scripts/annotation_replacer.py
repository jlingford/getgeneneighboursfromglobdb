#!/usr/bin/env python3

import pandas as pd

annotation_table = "../data/AMXMAG_0088_annotations_decompressed.tsv"
gff_input = "../data/AMXMAG_0088___3888___genetic_neighbourhood.gff"
new_gff_output = "../data/annotation_replacement_test2.gff"

# Load annotation table into a dict (indexed by ID, lightweight in memory)
ann_df = pd.read_csv(
    annotation_table,
    sep="\t",
    header=None,
    names=["ID", "db_xref", "Name", "product", "evalue"],
)

# Keep only rows where db_xref == "KOfam"
ann_df = ann_df[ann_df["db_xref"] == "KOfam"]
# print(ann_df.sort_values(by="ID"))

ann_df = ann_df.loc[ann_df.groupby("ID")["evalue"].idxmax()]
# print(ann_df.sort_values(by="ID"))

ann_map = ann_df.set_index("ID")[["db_xref", "Name", "product"]].to_dict(orient="index")

with open(gff_input) as infile, open(new_gff_output, "w") as outfile:
    for line in infile:
        if line.startswith("#"):
            # keep GFF header/comment lines intact
            outfile.write(line)
            continue

        parts = line.strip().split("\t")
        if len(parts) < 9:
            outfile.write(line)
            continue

        seqid, source, type_, start, end, score, strand, phase, attributes = parts

        # parse attributes into dict
        attr_dict = dict(
            field.split("=", 1) for field in attributes.split(";") if "=" in field
        )
        gene_id = attr_dict.get("ID")

        if gene_id and gene_id in ann_map:
            ann = ann_map[gene_id]
            # replace target fields
            attr_dict.update(
                {
                    "Name": ann["Name"],
                    "db_xref": ann["db_xref"],
                    "product": ann["product"],
                }
            )

        # rebuild attributes string
        new_attr = ";".join(f"{k}={v}" for k, v in attr_dict.items() if v != "")

        # write back GFF line
        outfile.write(
            "\t".join(
                [seqid, source, type_, start, end, score, strand, phase, new_attr]
            )
            + "\n"
        )
