# Gene neighbourhood extractor from GlobDB

This python script is designed to retrieve other genes surrounding a target gene of interest that are available within the GlobDB.

**About the GlobDB**: The GlobDB is a collection of species dereplicated genomes, combining and consolidating genomes from multiple different databases, such as the GTDB.
The GlobDB has a few big advantages over other publicly available genome collections.
First, every gene fasta file has a simple and straightforward naming system (e.g. GENOMEID___NUMBER), unlike the GTDB.
Second, it comes with useful accessory files such as HMM annotations for every gene and gff3 files. Both are very useful for plotting and analysing gene neighbourhoods.
