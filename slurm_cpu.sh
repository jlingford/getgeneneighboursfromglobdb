#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J james
#SBATCH --mem=40000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --account=rp24
#SBATCH --partition=genomicsb
#SBATCH --qos=genomicsbq
#SBATCH --time=48:00:00
#SBATCH --mail-user=james.lingford@monash.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_OUT
#SBATCH --error=log-%j.err
#SBATCH --output=log-%j.out

conda activate XXX

globdb_dir=
outputdir=
inputlist=
hyddbclasses=
arialfont=

python nife_ssu_finder.py \
    -d "$globdb_dir" \
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
    --svg_format
