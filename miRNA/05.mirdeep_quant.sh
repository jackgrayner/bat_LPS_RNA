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

cat novel_pha_rnas_filtered_precursor.fa all_hairpin.fa > novel_all_hairpin.fa
cat novel_pha_rnas_filtered_mature.fa mirbase_bats_mammal_mature.fa > novel_all_mature.fa

#remove duplicate mature sequences 
#switch to columns
awk 'NR%2==1 {h=$0} NR%2==0 {print $0"\t"h}' novel_all_mature.fa \
| sort -k1,1 \
| awk '{
    if (!seen[$1]) {
        name=$2
        sub(/^>/,"",name)
        print ">"name"\n"$1
        seen[$1]=1
    }
}' >  novel_all_mature_dedupe.fa

awk 'NR%2==1 {h=$0} NR%2==0 {print $0"\t"h}' novel_all_hairpin.fa \
| sort -k1,1 \
| awk '{
    if (!seen[$1]) {
        name=$2
        sub(/^>/,"",name)
        print ">"name"\n"$1
        seen[$1]=1
    }
}' >  novel_all_hairpin_dedupe.fa

#now run quantification
while read sample
do
	/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/quantifier.pl -y ${sample} -d -p novel_all_hairpin_dedupe.fa -m novel_all_mature_dedupe.fa -r ../${sample}_collapsed.fa 
done < ../../sample_list

#combine results files
while read sample
do
	cat miRNAs_expressed_all_samples_${sample}.csv | awk '{print $2}'  | sed "s/read_count/${sample}/g" > ${sample}_expression
done < ../../sample_list

cat miRNAs_expressed_all_samples_2311_TH.csv | awk '{print $1}' > mirna_ids
paste mirna_ids *_expression > all_samples_expression.tsv

