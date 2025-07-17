import pandas as pd
import dna_features_viewer as dfv
from bokeh.plotting import show, output_notebook
from bokeh.models import Range1d

# --- Configuration ---
gff_file = "nife_test.gff"
gene_name = "AMXMAG_0088___3611"  # <-- Replace with your gene name
window = 10000  # Window size upstream/downstream

output_notebook()

# --- Step 1: Locate gene coordinates from GFF ---
gff = pd.read_csv(
    gff_file,
    sep="\t",
    comment="#",
    header=None,
    names=[
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

# Find the gene by Name= or ID= (adjust as needed)
gene_row = gff[
    gff["attributes"].str.contains(f"Name={gene_name}|ID={gene_name}", na=False)
]

if gene_row.empty:
    raise ValueError(f"Gene '{gene_name}' not found in {gff_file}")


gene_start = gene_row.iloc[0]["start"]
gene_end = gene_row.iloc[0]["end"]
crop_start = gene_start - window
crop_end = gene_end + window


# --- Step 2: Translate to GraphicRecord ---
record = dfv.BiopythonTranslator().translate_record(gff_file)
print(record)

# --- Step 3: Manually crop and adjust features ---
cropped_features = []
for f in record.features:
    if f.end >= crop_start and f.start <= crop_end:
        new_start = max(f.start, crop_start) - crop_start
        new_end = min(f.end, crop_end) - crop_start
        new_feature = dfv.GraphicFeature(
            start=new_start,
            end=new_end,
            strand=f.strand,
            color=f.color,
            label=f.label,
            label_link_color=f.label_link_color,
            fontdict=f.fontdict,
            linewidth=f.linewidth,
            linecolor=f.linecolor,
        )
        cropped_features.append(new_feature)

# --- Step 4: Create cropped record and plot ---
cropped_record = dfv.GraphicRecord(
    sequence_length=crop_end - crop_start, features=cropped_features
)

fig = cropped_record.plot_with_bokeh()

# --- Step 5: Adjust Bokeh axis to show original genomic coordinates ---
fig.x_range = Range1d(start=crop_start, end=crop_end)
fig.xaxis.axis_label = "Genomic Position (bp)"
fig.xaxis[0].formatter.use_scientific = False  # Optional: avoid sci notation

# Optional: re-adjust label positions if they overlap due to new scale

# --- Step 6: Show the plot ---
show(fig)
