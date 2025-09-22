#!/usr/bin/env python3

import matplotlib.pyplot as plt
from dna_features_viewer import (
    GraphicFeature,
    GraphicRecord,
    BiopythonTranslator,
)


features = [
    GraphicFeature(start=5, end=20, strand=+1, color="#ffd700", label="Small feature"),
    GraphicFeature(
        start=20,
        end=500,
        strand=+1,
        color="#ffcccc",
        label="Gene 1 with a very long name",
    ),
    GraphicFeature(start=400, end=700, strand=-1, color="#cffccc", label="Gene 2"),
    GraphicFeature(start=600, end=900, strand=+1, color="#ccccff", label="Gene 3"),
]


# PLOT AND EXPORT A LINEAR VIEW OF THE CONSTRUCT
record = GraphicRecord(
    sequence_length=1000,
    features=features,
)
ax, _ = record.plot(
    figure_width=15,
    strand_in_label_threshold=3,
    annotate_inline=False,
    max_label_length=20,
    max_line_length=5,
    level_offset=1,
    draw_line=True,
    elevate_outline_annotations=True,
)
ax.figure.savefig("plot_linear.png", bbox_inches="tight")


# # PLOT AND EXPORT A CIRCULAR VIEW OF THE CONSTRUCT
# circular_rec = CircularGraphicRecord(sequence_length=1000, features=features)
# ax2, _ = circular_rec.plot(figure_width=4)
# ax2.figure.tight_layout()
# ax2.figure.savefig("plot_circular.png", bbox_inches="tight")


record = GraphicRecord(
    sequence=250 * "ATGC",
    features=[
        GraphicFeature(
            start=10, end=20, strand=+1, color="#ffd700", label="Small feature"
        ),
        GraphicFeature(
            start=20,
            end=500,
            strand=+1,
            color="#ffcccc",
            label="Gene 1 with a very long name",
        ),
        GraphicFeature(start=400, end=700, strand=-1, color="#cffccc", label="Gene 2"),
        GraphicFeature(start=600, end=900, strand=+1, color="#ccccff", label="Gene 3"),
    ],
)
zoom_start, zoom_end = 398, 428  # coordinates of the "detail"
cropped_record = record.crop((zoom_start, zoom_end))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 2))

# PLOT THE WHOLE SEQUENCE

ax1.set_title("Whole sequence", loc="left", weight="bold")
record.plot(ax=ax1)
ax1.fill_between((zoom_start, zoom_end), +1000, -1000, alpha=0.15)

# PLOT THE SEQUENCE DETAILS

cropped_record.plot_translation(
    ax=ax2, location=(408, 423), fontdict={"weight": "bold"}
)
cropped_record.plot(ax=ax2, plot_sequence=True)
ax2.set_title("Sequence detail", loc="left", weight="bold")

fig.savefig("plot_detail.png", bbox_inches="tight")
