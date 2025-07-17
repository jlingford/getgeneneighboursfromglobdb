#!/usr/bin/env python3

from dna_features_viewer import BiopythonTranslator
from dna_features_viewer import GraphicRecord


graphic_record = BiopythonTranslator().translate_record("./nife_test.gff")

# graphic_record = graphic_record.GraphicRecord()
graphic_record = graphic_record.crop((200000, 220000))

ax, _ = graphic_record.plot(figure_width=20, strand_in_label_threshold=3)
ax.figure.tight_layout()
ax.figure.savefig("test2.png")
