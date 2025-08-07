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
cat ../*collapsed.fa > all_collapsed.fa
/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/mapper.pl all_collapsed.fa -e -o 8 -h -i -j -m -s all_collapsed.fa -t all_collapsed.arf -p ../phha
/mnt/apps/users/jrayner/conda/envs/mirdeep2/bin/miRDeep2.pl all_collapsed.fa ../phha_genomic_nowhitespace.fa all_collapsed.arf none ../bats_mammal_mature_nowhitespace.fa ../mammal_hairpin_nowhitespace.fa 
