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

#remove duplicate mature sequences (sorted by seq, then sorted by name, then dupes removed - in R)
while read line
do
  grep -A 1 -w $line phha_mammal_mature.fa >> phha_mammal_mature_unique.fa
done < keep_dedup.txt

while read sample
do
	/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/quantifier.pl -y ${sample} -d -p ../phha_mammal_hairpin.fa -m ../phha_mammal_mature_unique.fa -r ../${sample}_collapsed.fa 
done < ../../sample_list

while read sample
do
	cat miRNAs_expressed_all_samples_${sample}.csv | awk '{print $2}' > ${sample}_expression
done < ../../sample_list

cat miRNAs_expressed_all_samples_2311_TH.csv | awk '{print $3}' > mirna_ids
paste mirna_ids *_expression > all_samples_expression.tsv

