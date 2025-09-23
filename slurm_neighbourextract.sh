#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J geneneighbours
#SBATCH --mem=40000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --account=rp24
#SBATCH --partition=genomics
#SBATCH --qos=genomics
#SBATCH --time=4:00:00
#SBATCH --mail-user=james.lingford@monash.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_OUT
#SBATCH --error=%j.err
#SBATCH --output=%j.out

conda activate /fs04/scratch2/rp24/jamesl2/MMseqs2/getgeneneighboursfromglobdb/rp24_scratch2/jamesl2/miniconda/conda/envs/geneneighbours

# take input args
inputdir=$1
count=$2

# set inputs for python script
# globdb_dir='/home/jamesl/rp24_scratch/Database/GlobDB_r226/globdb_r226'
outputdir="./output/nife_chunk${count}"
# inputlist='./input/nife_mmseqs_hits_IDs.txt'
inputlist="${inputdir}"
hyddbclasses='./auxfiles/nife_geneids2hyddb.tsv'
arialfont='./auxfiles/arial.ttf'
pickleindex='./globdb_index.pkl'

mkdir -p "$outputdir"

python nife_ssu_finder.py \
    -i "$pickleindex" \
    -l "$inputlist" \
    -o "$outputdir" \
    --add_hyddb "$hyddbclasses" \
    --add_arial_font "$arialfont" \
    --pairwise_set1 \
    --pairwise_set2 \
    --boltz_fastas \
    --colabfold_fastas \
    --chai_fastas \
    --dpi 100 \
    --pdf_format \
    --svg_format >output_${count}.log
