#!/bin/bash
#SBATCH --job-name=Kal_idx
#SBATCH -t 4:0:0
#SBATCH -n 16
#SBATCH --mem-per-cpu=8024
#SBATCH --oversubscribe
#SBATCH --output=Kal_idx.log
module load kallisto
spp=$1
kallisto index -i $spp *_rna.fna.gz
