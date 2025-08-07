#!/bin/bash
#SBATCH --job-name=mirdeep_map
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jrayner1@umd.edu
#SBATCH --output=mirdeep_map.log
source activate mirdeep2
while read sample
do
	/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/mapper.pl ../filtered_reads/${sample}_filtered.fastq -e -o 4 -h -i -j -l 18 -m -s ${sample}_collapsed.fa -t ${sample}_vs_genome.arf -p phha -v
done < ../sample_list
