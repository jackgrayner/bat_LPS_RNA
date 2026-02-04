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


setwd("~/Documents/UMD_new/nih_rna/baboon_reanalysis/")
dir<-"~/Documents/UMD_new/nih_rna/baboon_reanalysis/"
samples <- read.csv("SraRunTable.csv", header = TRUE)
# files <- file.path(paste0(dir,"/kallisto/", samples$Run, "/abundance.h5"))
# summary(all(file.exists(files)))
# 
# gtf <- file.path("GCF_008728515.1_Panubis1.0_genomic.gff.gz")
# txdb<-makeTxDbFromGFF(gtf,format="gff")
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
# txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
# cts<-txi$counts
# colnames(cts)<-samples$Run
# write.csv(cts,file="baboons_counts.csv")

cts<-read.csv("baboons_counts.csv",row.names = 1)

rownames(samples)<-samples$Run

cts<-cts[rowSums(cpm(cts)>0)>=nrow(samples)/4, ] #keep only genes expressed in > 1 in 4 samples
cts<-cts[,colnames(cts) %in% rownames(samples)] #
summary(rownames(samples)==colnames(cts)) #ensure sample order is the same

samples$ID<-gsub("_LPS","",samples$Sample.Name)
samples$ID<-gsub("_NULL","",samples$ID)
samples$ID<-as.factor(samples$ID)

dds.bab <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts),#remove unpaired samples, T and Z
  colData = samples,
  design= ~ treatment+ID))
resultsNames(dds.bab)

res.th.vs.c0<-add_gene_info(data.frame(results(dds.bab, name="treatment_NULL_vs_LPS",alpha=0.05)))
write.csv(res.th.vs.c0,"baboon_LPS.csv",row.names=FALSE,quote=FALSE)

bab.vsd<-assay(vst(dds.bab,blind=FALSE))
write.csv(vsd,"~/Documents/UMD_new/nih_rna/baboon_reanalysis/baboon_vst.csv",row.names=TRUE,quote=FALSE)
#vsd<-read.csv("~/Documents/UMD_new/nih_rna/baboon_reanalysis/baboon_vst.csv")

