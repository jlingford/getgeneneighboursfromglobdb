#!/usr/bin/env python3

from numpy import isin

import dna_features_viewer as dfv
from BCBio import GFF
# from dna_features_viewer import BiopythonTranslator
# from dna_features_viewer import GraphicRecord

# file = GFF.parse("./data/nife_test.gff")
file = "./data/nife_test.gff"

# graphic_record = dfv.BiopythonTranslator().translate_record("./data/nife_test.gff")
graphic_record = dfv.BiopythonTranslator().translate_record(file)

# graphic_record = graphic_record.GraphicRecord()
graphic_record = graphic_record.crop((200000, 220000))

ax, _ = graphic_record.plot(figure_width=20, strand_in_label_threshold=3)
ax.figure.tight_layout()
ax.figure.savefig("output/test2.png")

#################################################################

from dna_features_viewer import BiopythonTranslator
import matplotlib.pyplot as plt

# Your target gene of interest
gff_input_file = "./data/AMXMAG_0088___3888___genetic_neighbourhood.gff"
target_gene_id = "3888"
# window_start = 200000
# window_end = 220000
# anno_code_list = ["COG1008", "COG1009"]
anno_codes_neighbours = [
    "COG1941",
    "COG1908",
    "COG5557",
    "COG1740",
]  # These are the FrhA, FrhG, FrhD subunits respectively for NiFe hydrogenase
anno_codes_targets = [
    "COG0374",
    "COG3259",
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
                if gene_anno in anno_codes_neighbours:
                    return "#89b4fa"
                # colour expected target with correct annotation as green
                if gene_anno in anno_codes_targets:
                    return "#a6e3a1"
            # colour target gene red if it doesn't match anything else
            if target_gene_id == gene_id:
                return "#f38ba8"
        return "#cdd6f4"


# class CustomTranslator(BiopythonTranslator):
#     def compute_feature_color(self, feature):
#         # Color the target gene red, all others light grey
#         if "ID" in feature.qualifiers:
#             print(feature.qualifiers)
#             gene_id = feature.qualifiers["ID"][0]
#             if "Name" in feature.qualifiers:
#                 gene_anno = feature.qualifiers["Name"][0]
#                 if gene_anno in anno_code_list:
#                     return "#89b4fa"  # Target gene color
#                 print(gene_anno)
#             print(gene_id)
#             gene_id = gene_id.split("___")[1]
#             # if gene_id == target_gene_id:
#             if target_gene_id == gene_id:
#                 return "#a6e3a1"  # Target gene color
#         return "#cdd6f4"  # Default color for other genes


# Translate using the custom translator
record = CustomTranslator().translate_record(gff_input_file)

# Set sequence length to avoid cropping errors
# record.sequence_length = max(window_end, record.sequence_length)

# Crop to window region
# record = record.crop((window_start, window_end))

# Plot figure

fig, ax = plt.subplots()
ax, _ = record.plot(figure_width=20, strand_in_label_threshold=0)
ax.figure.tight_layout()
ax.figure.savefig("output/color_test.png")
plt.close("all")

###########################################################


# def compute_feature_color(feature):
#     # Color the target gene red, all others light grey
#     if "ID" in feature.qualifiers:
#         print(feature.qualifiers)
#         gene_id = feature.qualifiers["ID"][0]
#         print(gene_id)
#         # if gene_id == target_gene_id:
#         if target_gene_id in gene_id:
#             print(target_gene_id)
#             return "#a6e3a1"  # Target gene color
#     return "#a6adc8"  # Default color for other genes
#
#
# compute_feature_color(gff_input_file)
