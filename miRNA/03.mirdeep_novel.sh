#!/bin/bash
#SBATCH --job-name=mirdeep_novel
#SBATCH --export=ALL
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=long
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jrayner1@umd.edu
#SBATCH --output=mirdeep_novel.log
source activate mirdeep2
#concatenate
cat ../*collapsed.fa > all_collapsed1.fa
#collapse concatenated fasta
/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/mapper.pl all_collapsed1.fa -c -o 8 -i -j -m -s all_collapsed.fa -t all_collapsed.arf -p ../phha

#keep only human, mouse, bat miRNAs
grep -A 1 "^>hsa\|^>mmu\|^>aja\|^>pal\|^>efu" hairpin_nowhitespace.fa | grep -v "\-\-" > mirbase_bats_mammal_hairpin.fa
grep -A 1 "^>hsa\|^>mmu\|^>aja\|^>pal\|^>efu" mature_nowhitespace.fa | grep -v "\-\-" > mirbase_bats_mammal_mature.fa
cat mirbase_bats_mammal_hairpin.fa mlu_precursor.fa mmy_precursor.fa > all_hairpin.fa

#now identify novel transcripts
/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/miRDeep2.pl all_collapsed.fa ../phha_genomic_nowhitespace.fa all_collapsed.arf none mirbase_bats_mammal_mature.fa all_hairpin.fa
