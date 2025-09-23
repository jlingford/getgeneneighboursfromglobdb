#!/bin/bash -l
#SBATCH -D ./
#SBATCH -J james
#SBATCH --mem=10000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --account=rp24
#SBATCH --partition=genomics
#SBATCH --qos=genomics
#SBATCH --time=4:00:00
#SBATCH --mail-user=james.lingford@monash.edu
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_OUT
#SBATCH --error=%j.err
#SBATCH --output=%j.out

# ---
INPUTDIR='./input/split_input_lists'

count=-1
for file in ${INPUTDIR}/*; do
    if [[ -f "$file" ]]; then
        ((count++))
        n=$(printf %03d $count)
        echo $n
        sbatch -J nife_ssu_${n} slurm_neighbourextract.sh "$file" $n
    fi
done
