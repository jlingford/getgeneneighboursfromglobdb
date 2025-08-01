import tarfile
import os
import argparse

def load_genome_ids(genome_list_file):
    with open(genome_list_file) as f:
        return set(line.strip() for line in f if line.strip())

def match_target(filename, genome_ids):
    base = os.path.basename(filename)
    for gid in genome_ids:
        if base == f"{gid}.fa.gz":
            return True
    return False

def extract_matching_files(tar_path, genome_ids, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with tarfile.open(tar_path, "r:gz") as tar:
        for member in tar.getmembers():
            if not member.isfile():
                continue
            if match_target(member.name, genome_ids):
                print(f"Extracting: {member.name}")
                extracted = tar.extractfile(member)
                out_path = os.path.join(output_dir, os.path.basename(member.name))
                with open(out_path, 'wb') as out_f:
                    out_f.write(extracted.read())

def main():
    parser = argparse.ArgumentParser(description="Extract specific .fa.gz files from a tar.gz archive.")
    parser.add_argument("archive", help="Path to the .tar.gz archive")
    parser.add_argument("genome_list", help="Text file with genome IDs (one per line)")
    parser.add_argument("-o", "--output", default="extracted_fastas", help="Output directory (default: extracted_fastas)")
    args = parser.parse_args()

    genome_ids = load_genome_ids(args.genome_list)
    extract_matching_files(args.archive, genome_ids, args.output)

if __name__ == "__main__":
    main()
