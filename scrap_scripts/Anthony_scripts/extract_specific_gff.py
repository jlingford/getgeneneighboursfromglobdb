import tarfile
import gzip
import os
import argparse

def load_genome_ids(genome_list_file):
    with open(genome_list_file) as f:
        return set(line.strip() for line in f if line.strip())

def match_target(filename, genome_ids):
    for gid in genome_ids:
        if filename.endswith(f"{gid}_cog.gff.gz"):
            return gid
    return None

def extract_matching_files(tar_path, genome_ids, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with tarfile.open(tar_path, "r:gz") as tar:
        members = tar.getmembers()
        for member in members:
            if not member.isfile():
                continue
            matched_id = match_target(os.path.basename(member.name), genome_ids)
            if matched_id:
                print(f"Extracting: {member.name}")
                file_content = tar.extractfile(member)
                out_path = os.path.join(output_dir, os.path.basename(member.name))
                with open(out_path, 'wb') as out_f:
                    out_f.write(file_content.read())

def main():
    parser = argparse.ArgumentParser(description="Extract specific GFF files from a tar.gz archive.")
    parser.add_argument("archive", help="Path to .tar.gz file")
    parser.add_argument("genome_list", help="Path to text file with genome IDs (one per line)")
    parser.add_argument("-o", "--output", default="extracted_gffs", help="Output directory (default: extracted_gffs)")
    args = parser.parse_args()

    genome_ids = load_genome_ids(args.genome_list)
    extract_matching_files(args.archive, genome_ids, args.output)

if __name__ == "__main__":
    main()
