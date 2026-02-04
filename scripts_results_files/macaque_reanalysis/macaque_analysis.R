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


setwd("~/Documents/UMD_new/nih_rna/macaque_reanalysis")
# dir<-"~/Documents/UMD_new/nih_rna/macaque_reanalysis"
# param<-"180_20"
# param<-"200_30"
samples <- read.csv("SraRunTable.csv", header = TRUE)
# files <- file.path(paste0(dir,"kallisto_", param,"/",samples$Run, "/abundance.h5"))
# summary(all(file.exists(files)))
# 
# gtf <- file.path("GCF_003339765.1_Mmul_10_genomic.gff.gz")
# txdb<-makeTxDbFromGFF(gtf,format="gff")
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
# txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
# cts<-txi$counts
# colnames(cts)<-samples$Run
# write.csv(cts,file=paste0("macaques_counts_",param,".csv"))

cts<-read.csv("macaques_counts_180_20.csv",h=T,row.names = 1)
rownames(samples)<-samples$Run
cts<-cts[rowSums(cpm(cts)>0)>=round(nrow(samples)/4), ] #keep only genes expressed in > 1 in 4 samples
cts<-cts[,colnames(cts) %in% rownames(samples)] #
summary(rownames(samples)==colnames(cts)) #ensure sample order is the same
samples$subject_id<-factor(samples$subject_id)
dds.mac <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts),#remove unpaired samples, T and Z
  colData = samples,
  design= ~ treatment+subject_id))
resultsNames(dds.mac)
res.th.vs.c0.1<-add_gene_info(data.frame(results(dds.mac, name="treatment_NC_vs_LPS",alpha=0.05)))

# cts<-read.csv("macaques_counts_200_30.csv",h=T,row.names = 1)
# rownames(samples)<-samples$Run
# cts<-cts[rowSums(cpm(cts)>0)>=round(nrow(samples)/4), ] #keep only genes expressed in > 1 in 4 samples
# cts<-cts[,colnames(cts) %in% rownames(samples)] #
# summary(rownames(samples)==colnames(cts)) #ensure sample order is the same
# samples$subject_id<-factor(samples$subject_id)
# dds <- DESeq(DESeqDataSetFromMatrix(
#   countData = round(cts),#remove unpaired samples, T and Z
#   colData = samples,
#   design= ~ treatment+subject_id))
# resultsNames(dds)
# res.th.vs.c0.2<-add_gene_info(data.frame(results(dds, name="treatment_NC_vs_LPS",alpha=0.05)))

# both<-merge(res.th.vs.c0.1,res.th.vs.c0.2,by='gene')
# cor.test(both$log2FoldChange.x,both$log2FoldChange.y)
# plot(both$log2FoldChange.x,both$log2FoldChange.y)
# 
write.csv(res.th.vs.c0.1,"macaque_LPS.csv",row.names=FALSE,quote=FALSE)


mac.vsd<-assay(vst(dds.mac,blind=FALSE))
write.csv(vsd,"~/Documents/UMD_new/nih_rna/macaque_reanalysis/macaque_vst.csv",row.names=TRUE,quote=FALSE)
