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


class CustomTranslator(BiopythonTranslator):
    def compute_feature_color(self, feature):
        # Color the target gene red, all others light grey
        if "ID" in feature.qualifiers:
            print(feature.qualifiers)
            gene_id = feature.qualifiers["ID"][0]
            if "Name" in feature.qualifiers:
                gene_anno = feature.qualifiers["Name"][0]
                if gene_anno in anno_code_list:
                    return "#89b4fa"  # Target gene color
                print(gene_anno)
            print(gene_id)
            gene_id = gene_id.split("___")[1]
            # if gene_id == target_gene_id:
            if target_gene_id == gene_id:
                return "#a6e3a1"  # Target gene color
        return "#cdd6f4"  # Default color for other genes


# Your target gene of interest
gff_input_file = "./data/nife_test.gff"
target_gene_id = "3611"
window_start = 200000
window_end = 220000
anno_code_list = ["COG1008", "COG1009"]

# Translate using the custom translator
record = CustomTranslator().translate_record(gff_input_file)

# Set sequence length to avoid cropping errors
record.sequence_length = max(window_end, record.sequence_length)

# Crop to window region
record = record.crop((window_start, window_end))

# Plot figure
ax, _ = record.plot(figure_width=20, strand_in_label_threshold=0)
ax.figure.tight_layout()
ax.figure.savefig("output/color_test.png")

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
