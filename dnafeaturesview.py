#!/usr/bin/env python3

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
