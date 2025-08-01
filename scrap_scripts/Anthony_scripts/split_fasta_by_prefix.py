from Bio import SeqIO
import sys
import os

def split_fasta_by_prefix(input_fasta):
    # Dictionary to keep file handles open per genome ID
    file_handles = {}

    # Iterate over each record in the input fasta
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.id  # This is the header without '>'
        if "___" not in header:
            print(f"Warning: skipping record with unexpected header format: {header}")
            continue
        genome_id = header.split("___")[0]

        # Open file handle if not already open
        if genome_id not in file_handles:
            out_filename = f"{genome_id}.faa"
            file_handles[genome_id] = open(out_filename, "w")

        # Write the record to the appropriate file
        SeqIO.write(record, file_handles[genome_id], "fasta")

    # Close all file handles
    for fh in file_handles.values():
        fh.close()

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python split_fasta_by_prefix.py input.faa")
        sys.exit(1)

    input_fasta = sys.argv[1]
    if not os.path.exists(input_fasta):
        print(f"Error: File '{input_fasta}' not found.")
        sys.exit(1)

    split_fasta_by_prefix(input_fasta)
    print("Splitting complete.")
