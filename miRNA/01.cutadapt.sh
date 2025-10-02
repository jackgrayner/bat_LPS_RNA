#!/bin/bash
#SBATCH --job-name=cutadapt
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jrayner1@umd.edu
source activate cutadapt
while read sample
do
	cutadapt -m 18 -M 26 -o ./${sample}_filtered.fastq.gz ../${sample}/ILLUMINA_DATA/*R1_trimmed.fastq.gz
	gunzip ${sample}_filtered.fastq.gz
done < ../sample_list
