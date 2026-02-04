library(viridis)
library(DESeq2) 
library(edgeR) 
library(pheatmap)
library(RColorBrewer)
library(ggplot2) 
library(reshape2) 
library(tximportData) #only needed on first run
library(tximport)
library(GenomicFeatures)
library(patchwork)
library(dplyr) 
library(clusterProfiler)
organism<-'org.Hs.eg.db'
library(organism, character.only = TRUE)

#define function to add extra gene info to results files
add_gene_info<-function(x){
  #x<-x[!is.na(x$padj),]#get rid of NAs
  x$sig<-FALSE
  x[x$padj<0.05 & !is.na(x$padj),]$sig<-TRUE
  x$gene<-rownames(x)
  return(x)
}


setwd("~/Documents/UMD_new/nih_rna/pigs/")
dir<-"~/Documents/UMD_new/nih_rna/pigs/"
samples <- read.csv("sample_info.csv", header = TRUE)
files <- file.path(paste0(dir,"/kallisto/", samples$Sample, "/abundance.h5"))
summary(all(file.exists(files)))

# gtf <- file.path("GCF_000003025.6_Sscrofa11.1_genomic.gff.gz")
# txdb<-makeTxDbFromGFF(gtf,format="gff")
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
# txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
# cts<-txi$counts
# colnames(cts)<-samples$Sample
# write.csv(cts,file="pigs_counts.csv")
cts<-read.csv("pigs_counts.csv",row.names = 1)

rownames(samples)<-samples$Sample

cts<-cts[rowSums(cpm(cts)>0)>=nrow(samples)/4, ] #keep only genes expressed in > 1 in 4 samples
cts<-cts[,colnames(cts) %in% rownames(samples)] #
summary(rownames(samples)==colnames(cts)) #ensure sample order is the same

samples$ID<-as.factor(samples$Characteristics.individual.)

dds.pig <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts),#remove unpaired samples, T and Z
  colData = samples,
  design= ~ treatment+ID))
resultsNames(dds)

res.th.vs.c0<-add_gene_info(data.frame(results(dds.pig, name="treatment_LPS_vs_C",alpha=0.05)))
res.th.vs.c0<-res.th.vs.c0[order(res.th.vs.c0$pvalue),]
write.csv(res.th.vs.c0,"pig_LPS.csv",row.names=FALSE,quote=FALSE)

pig.vsd<-assay(vst(dds.pig,blind=FALSE))
write.csv(vsd,"pig_vst.csv",row.names=TRUE,quote=FALSE)
