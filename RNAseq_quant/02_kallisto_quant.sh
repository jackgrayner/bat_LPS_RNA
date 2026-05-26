#!/bin/bash
#SBATCH --job-name=Kal_aln
#SBATCH -t 1:0:0
#SBATCH -n 1
#SBATCH --mem-per-cpu=8G
#SBATCH --oversubscribe
#SBATCH --array=1-28
#SBATCH --output=Kallisto_aln.%A_%a.out

#load the kallisto module
module load kallisto

#read spp name from input argument
spp=$1

#quantify samples S1 through S28 (or through arrange range defined above)
kallisto quant -i $spp -b 50 -o ./kallisto/S${SLURM_ARRAY_TASK_ID} ./*S${SLURM_ARRAY_TASK_ID}_R1*.fastq.gz ./*S${SLURM_ARRAY_TASK_ID}_R2*.fastq.gz 
