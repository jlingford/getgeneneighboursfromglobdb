#!/usr/bin/env python3

from Bio import SeqIO
from BCBio import GFF
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Convert .fa + .gff + .faa to GenBank (.gbk)"
    )
    parser.add_argument(
        "-f", "--fasta", required=True, help="Genome nucleotide fasta file (.fa)"
    )
    parser.add_argument("-g", "--gff", required=True, help="GFF annotation file (.gff)")
    parser.add_argument(
        "-a", "--faa", required=True, help="Translated protein file (.faa)"
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output GenBank file (.gbk)"
    )
    return parser.parse_args()


def load_proteins(faa_file):
    protein_seqs = {}
    for record in SeqIO.parse(faa_file, "fasta"):
        protein_seqs[record.id] = str(record.seq)
    return protein_seqs


def extract_gff_id(feature):
    if "ID" in feature.qualifiers:
        return feature.qualifiers["ID"][0]
    for key in ["id", "Name"]:
        if key in feature.qualifiers:
            return feature.qualifiers[key][0]
    return None


def flatten_features(rec):
    """
    Flatten all CDS features to top-level features in rec.features
    """
    new_features = []
    for feat in rec.features:
        # If feature has sub_features (like gene or mRNA), promote CDS subfeatures
        if hasattr(feat, "sub_features") and feat.sub_features:
            for subf in feat.sub_features:
                if subf.type == "CDS":
                    new_features.append(subf)
                else:
                    # Keep non-CDS subfeatures if you want, or skip
                    pass
        elif feat.type == "CDS":
            new_features.append(feat)
        else:
            # Keep non-CDS top-level features as well (e.g., genes)
            new_features.append(feat)
    rec.features = new_features


def gff_to_genbank(gff_file, fasta_file, faa_file, output_file):
    protein_seqs = load_proteins(faa_file)
    missing_translations = 0
    total_cds = 0

    with open(output_file, "w") as out_handle:
        with open(gff_file) as gff_handle, open(fasta_file) as fasta_handle:
            fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_handle, "fasta"))
            for rec in GFF.parse(gff_handle, base_dict=fasta_dict):
                rec.annotations["molecule_type"] = "DNA"
                flatten_features(rec)  # Make sure CDS are top-level

                for feat in rec.features:
                    if feat.type == "CDS":
                        total_cds += 1
                        qualifiers = feat.qualifiers
                        cds_id = extract_gff_id(feat)

                        if "product" not in qualifiers:
                            qualifiers["product"] = ["hypothetical protein"]

                        if cds_id:
                            qualifiers["locus_tag"] = [cds_id]
                            if cds_id in protein_seqs:
                                qualifiers["translation"] = [protein_seqs[cds_id]]
                            else:
                                missing_translations += 1
                        else:
                            missing_translations += 1

                SeqIO.write(rec, out_handle, "genbank")

    if total_cds == 0:
        print("⚠️ WARNING: No CDS features found in GFF file.")
    elif missing_translations > 0:
        print(f"⚠️ WARNING: {missing_translations} CDS features missing translations.")
    else:
        print(f"✅ All {total_cds} CDS translations included.")


if __name__ == "__main__":
    args = parse_arguments()
    gff_to_genbank(args.gff, args.fasta, args.faa, args.output)
