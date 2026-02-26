## miRNA analysis

- 01.cutadapt.sh: size filter miRNA reads
- 02.mirdeep_map.sh: map reads to P. hastatus genome
- 03.mirdeep_novel: identify novel miRNAs
- 04.filter_novel_miRNAs.R: filter novel miRNAs
- 05.mirdeep_quant.sh: quantify abundance of final miRNA set
- miRNA_analysis.R: perform differential abundance analysis in DESeq2
- novel_...: novel hairpin and mature sequences identified by miRDeep2
