#!/usr/bin/bash

annotation_table="../data/AMXMAG_0088_annotations_decompressed.tsv"
gff_input="../data/AMXMAG_0088___3888___genetic_neighbourhood.gff"
new_gff_output="../data/annotation_replacement_test.gff"

while IFS= read -r line; do
    gene_id=$(echo "$line" | awk -F"\t" '{print $9}' | awk -F";" '{print $1}' | sed 's/ID=//')
    partial_gff=$(echo "$line" | awk -F"\t" -vOFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8}')

    kofam_match=$(grep -P "${gene_id}\tKOfam" $annotation_table)

    db_xref=$(echo "$kofam_match" | awk -F"\t" '{print $2}')
    name=$(echo "$kofam_match" | awk -F"\t" '{print $3}')
    product=$(echo "$kofam_match" | awk -F"\t" '{print $4}')

    new_gff_attribute="ID=${gene_id};Name=${name};db_xref=${db_xref};product=${product}"
    new_gff_line="${partial_gff}\t${new_gff_attribute}"
    new_gff_line=$(echo "$new_gff_line" | sed 's/;Name=;db_xref=;product=//')

    echo -e "$new_gff_line"

done <$gff_input >$new_gff_output

# $new_gff_output
