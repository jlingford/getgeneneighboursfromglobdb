# fasta pairwise generator for Chai-1 input

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from itertools import combinations_with_replacement
from pathlib import Path


def generate_pairwise_fastas(
    input_fasta: Path, output_dir: Path, target_id: str
) -> None:
    """Generate multiple fasta files containing pairwise combinations of all gene gene_neighbours"""
    # params
    target_id = target_id
    # Read all sequences into a list
    records = list(SeqIO.parse(input_fasta, "fasta"))

    # Move the target gene of interest to the front
    # records.sort(key=lambda record: 0 if record.id == target_id else 1)

    # reorder all records by how adjacent they are to the GOI, with closest genes being listed first over more distant genes
    # get gene ids numbers for fasta ID first
    # gene_nums = [int(record.id.split("___")[1]) for record in records]
    # print(gene_nums)
    def dist_from_goi(record: SeqRecord):
        goi_num = int(target_id.split("___")[1])
        neighbour_num = int(record.id.split("___")[1])
        return abs(neighbour_num - goi_num)

    # sorted_records = sorted(records, key=lambda record: abs(gene_nums - int(target_id)) for  )
    # sorted_records = sorted(
    #     records,
    #     key=lambda record: abs(
    #         int(record.id.split("___")[1]) - int(target_id.split("___")[1])
    #     ),
    # )
    sorted_records = sorted(records, key=dist_from_goi)
    print(sorted_records)

    # Build a dictionary for quick access by ID (optional, if you want)
    id_to_record = {record.id: record for record in records}

    # Modify fasta header IDs for Chai-1 input
    renamed_records = [
        SeqRecord(
            seq=record.seq, name=record.id, id=f"protein|{record.id}", description=""
        )
        for record in sorted_records
    ]
    # print(renamed_records)

    # Generate all pairwise combinations (with replacement) using itertools
    pairs = combinations_with_replacement(renamed_records, 2)
    for i, (rec1, rec2) in enumerate(pairs, start=1):
        # Check if target is in this pair
        includes_target = target_id in {rec1.name, rec2.name}
        assert isinstance(rec1.name, str)
        assert isinstance(rec2.name, str)

        # get gene names to write new file names
        genome = rec1.name.split("___")[0]
        name1 = rec1.name.split("___")[1]
        name2 = rec2.name.split("___")[1]

        # Build output file name
        if includes_target:
            fname = f"{i:03d}_{genome}___{name1}-{name2}_geneofinterest.faa"
        else:
            fname = f"{i:03d}_{genome}___{name1}-{name2}.faa"

        print(fname)

        output_file = output_dir / fname

        # Write the pair
        SeqIO.write([rec1, rec2], output_file, "fasta")


# Example usage
input_fasta = Path(
    "./output/test2/AMXMAG_0088___3611/AMXMAG_0088___3611___gene_neighbours.faa"
)
gene_name = "AMXMAG_0088___3611"
output_dir = "./output/test3"
output_dir = Path(output_dir) / f"Fasta_pairs___{gene_name}"
output_dir.mkdir(parents=True, exist_ok=True)

generate_pairwise_fastas(input_fasta, output_dir, gene_name)

print(input_fasta.stem)

gene_name = input_fasta.stem.split("___")[1]
print(gene_name)
