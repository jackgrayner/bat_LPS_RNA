#!/bin/bash
#SBATCH --job-name=mirdeep_quant
#SBATCH --export=ALL
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jrayner1@umd.edu
#SBATCH --output=mirdeep_quant.log
source activate mirdeep2
while read sample
do
	/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/quantifier.pl -y ${sample} -d -p ../phha_mammal_hairpin.fa -m ../phha_mammal_mature.fa -r ../${sample}_collapsed.fa 
done < ../../sample_list
