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
cat ../phha_novel_mature.fa ../mammal_mature_nowhitespace.fa > phha_mammal_mature.fa
cat ../phha_novel_precusor.fa ../mammal_hairpin_nowhitespace.fa > phha_mammal_hairpin.fa

while read sample
do
	/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/quantifier.pl -p phha_mammal_hairpin.fa -m /phha_mammal_mature.fa -r ../${sample}_collapsed.fa 
done < ../../sample_list
