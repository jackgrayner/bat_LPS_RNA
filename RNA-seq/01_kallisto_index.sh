#!/bin/bash
#SBATCH --job-name=Kal_idx
#SBATCH -t 4:0:0
#SBATCH -n 16
#SBATCH --mem-per-cpu=8024
#SBATCH --oversubscribe
#SBATCH --output=Kal_idx.log

#load the kallisto module
module load kallisto

#read species name from input argument
spp=$1

#index the rna fasta file downloaded from NCBI
kallisto index -i $spp *_rna.fna.gz
