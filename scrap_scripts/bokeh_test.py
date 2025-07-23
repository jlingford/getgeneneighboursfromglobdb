import dna_features_viewer as dfv
import bokeh
import bokeh.plotting
# bokeh.io.output_notebook()

features = [
    dfv.GraphicFeature(
        start=5, end=20, strand=+1, color="#ffd700", label="Small feature"
    ),
    dfv.GraphicFeature(
        start=20,
        end=500,
        strand=+1,
        color="#ffcccc",
        label="Gene 1 with a very long name",
    ),
    dfv.GraphicFeature(start=400, end=700, strand=-1, color="#cffccc", label="Gene 2"),
    dfv.GraphicFeature(start=600, end=900, strand=+1, color="#ccccff", label="Gene 3"),
]
for f in features:
    f.linecolor = "green"
    f.box_color = "orange"  # box color is not being changed??
    f.label_link_color = "#ccccff"  # label link color is note being changed??

record = dfv.GraphicRecord(sequence_length=1000, features=features)
record.plot_with_bokeh()

# bokeh.plotting # bokeh figure is not being output to notebook because it is formatted as html doc
