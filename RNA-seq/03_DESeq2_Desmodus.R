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

# first run only - read in Kallisto files and sum transcripts across genes, save as CSV
dir<-paste0("~/Documents/UMD_new/desmodus_rna/")
samples <- read.csv("Dero_sample_info.csv", header = TRUE)
files <- file.path(paste0(dir, "kallisto/", samples$Sample, "/abundance.h5"))
summary(all(file.exists(files)))
gtf <- file.path("GCF_022682495.2_HLdesRot8A.1_genomic.gff.gz")
txdb<-makeTxDbFromGFF(gtf,format="gff")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")
txi <- tximport(files, type = "kallisto", tx2gene = tx2gene)
cts<-txi$counts
colnames(cts)<-samples$ID
write.csv(cts,file=paste0("cts_Dero",".csv"))

#####################
#FUNCTIONS AND THEMES
#####################

#theme for plots
SexPalette<-c("#C74955","#5090bf")
cust.theme<-function(){
  theme_minimal() %+replace%
    theme(panel.grid.minor=element_blank(),
          plot.background=element_rect(fill='white',colour='white'),
          panel.background=element_rect(fill="#fcfbfa",colour='black',linewidth = 0.25))
}

#define function to add extra gene info to DESeq results files
loc.ph.ids<-read.csv("LOC_ph_dr.csv",h=T) #read in reciprocal blastn best hits for ncRNAs
add_gene_info<-function(x){
  x$sig<-FALSE
  x[x$padj<0.05 & !is.na(x$padj),]$sig<-TRUE #nb. we retain NAs for overrepresentation tests
  x$gene<-rownames(x)
  x$note<-""
  x[x$gene %in% loc.ph.ids$LOC.dr,]$note<-"LOC ID replaced by PhHa LOC Id"
  y<-merge(x[x$gene %in% loc.ph.ids$LOC.dr,],loc.ph.ids,by.x="gene",by.y="LOC.dr")
  y$gene<-y$LOC.ph
  y<-y[,-10]
  rownames(y)<-y$gene
  length(unique(loc.ph.ids$LOC.dr))
  x<-x[!x$gene %in% loc.ph.ids$LOC.dr,]
  x<-rbind(x,y)
  return(x)
}

#define function to plot GO overrepresentation
plot_GO_overrep<-function(results_table,of_interest){
  go.OR<-enrichGO(gene = of_interest$gene,
                  universe = results_table$gene,#list of all genes
                  keyType = "SYMBOL",
                  OrgDb = organism,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  readable = TRUE)
  return(dotplot(go.OR)+cust.theme()+scale_fill_viridis(option="A",direction=-1)+
           theme(axis.ticks.y=element_line()))
}

plot_kegg_enrich<-function(results_table,genes){
  ids<-bitr(genes,fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)#convert symbol to ENTREZID (human)
  ids<-ids[!duplicated("ENTREZID"),]
  df2 = results_table[results_table$gene %in% ids$SYMBOL,]
  df2$SYMBOL<-df2$gene
  df2<-merge(df2,ids,by='SYMBOL')
  df2<-df2[order(df2$diff,decreasing = TRUE),]
  kegg_gene_list <- df2$diff
  names(kegg_gene_list) <- df2$ENTREZID
  geneList = sort(kegg_gene_list, decreasing = TRUE)
  kk2 <- gseKEGG(geneList     = (geneList),
                 organism     = 'human',
                 minGSSize    = 3,
                 maxGSSize    = 800,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",scoreType='pos',
                 keyType       = "ncbi-geneid")
  return(dotplot(kk2)+cust.theme()+scale_fill_viridis(direction=-1,option="A"))
}

#plot volcano
plot_volcano<-function(df){
  if(nrow(df[df$pvalue==0 & !is.na(df$pvalue),])>0){
    df[df$pvalue==0 & !is.na(df$pvalue),]$pvalue<-1.520048e-299
    }
  df<-df[!is.na(df$padj),]
  anno.set<-df[abs(df$log2FoldChange)>quantile(abs(df$log2FoldChange),0.95) | 
                 -log10(df$pvalue)>quantile(-log10(df$pvalue),0.95),]
  anno.set<-anno.set[-grep("^LOC",anno.set$gene),]
  g.volcano<-ggplot(df,aes(x=log2FoldChange,y=-log10(pvalue)))+cust.theme()+
    geom_vline(xintercept=0,linetype='dashed')+
    geom_point(aes(colour=padj<0.05),size=2,alpha=0.5)+scale_colour_manual(values=c("#aaaaaa","#e07575"))+
    geom_text_repel(data=anno.set,colour='black',fontface = 'italic',
                    aes(label=gene),size=2.5,max.overlaps = 8,min.segment.length = 0)
  return(g.volcano)
}

#############

setwd("~/Documents/UMD_new/desmodus_rna/")
samples <- read.csv("Dero_sample_info.csv", header = TRUE)
rownames(samples)<-samples$ID

cts<-read.csv("cts_Dero.csv",row.names=1)
colnames(cts)<-gsub("X","",colnames(cts))
cts<-cts[rowSums(cpm(cts)>0)>=nrow(samples)/4, ] #keep only genes expressed in > 1 in 4 samples
cts<-cts[,colnames(cts) %in% rownames(samples)] #
summary(rownames(samples)==colnames(cts)) #ensure sample order is the same

samples$Est.Age<-scale(as.numeric(samples$Age_may25)) #z-scale age variable for spp. comparability
samples$Band<-as.factor(samples$Band) #ensure band is treated as a factor 
dds <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts)[,!samples$Band %in% c("T","Z")],#remove unpaired samples, T and Z
  colData = samples[!samples$Band %in% c("T","Z"),],
  design= ~ Trt+Band))
resultsNames(dds)

#plot PCA of  1000 most-variable genes
pcaData <- plotPCA(vst(dds),intgroup=c("ID","Trt","Sex","Band","Est.Age"), 
                   returnData=TRUE,ntop=1000)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pca.all.vsd<-ggplot(pcaData, aes(x=PC1, y=PC2, colour=Trt)) +
  cust.theme()+
  geom_point(size=2)+geom_text_repel(aes(label=ID))+
  xlab(paste0("PC1 (",percentVar[1],"%)"))+
  ylab(paste0("PC2 (",percentVar[2],"%)"))

#plot PC1 diffs (separating treatments) vs age
pc1.comp<-merge(pcaData[pcaData$Trt=="C0",],
                pcaData[pcaData$Trt=="LPS",],by='Band')#create wide format df
pc1.comp$diff<-pc1.comp$PC1.x-pc1.comp$PC1.y
ggplot(pc1.comp[pc1.comp$Sex.x=="F",],aes(x=Est.Age.x,y=diff))+
  cust.theme()+geom_point()

#save DEG results
summary((results(dds, name="Trt_LPS_vs_C0",alpha=0.05)))
res.th.vs.c0<-add_gene_info(data.frame(results(dds, name="Trt_LPS_vs_C0",alpha=0.05)))
write.csv(res.th.vs.c0,"DeRo_TH_vs_C0.csv",row.names=FALSE,quote=FALSE)

#compare results across species
ph.th.vs.c0<-read.csv("../nih_rna/july_25/results_files/phase7_paired_th_vs_c0.csv",h=T)
dr.ph<-merge(res.th.vs.c0,ph.th.vs.c0,by="gene")
colnames(dr.ph)<-gsub("\\.x",".DeRo",colnames(dr.ph))
colnames(dr.ph)<-gsub("\\.y",".PhHa",colnames(dr.ph))
summary(factor(paste(dr.ph$sig.DeRo,dr.ph$sig.PhHa)))

PearsonsR<-cor.test(dr.ph$log2FoldChange.DeRo,dr.ph$log2FoldChange.PhHa)[,] #Pearson's correlation across genes
cor.test(dr.ph[dr.ph$sig.DeRo | dr.ph$sig.PhHa,]$log2FoldChange.DeRo,
         dr.ph[dr.ph$sig.DeRo | dr.ph$sig.PhHa,]$log2FoldChange.PhHa) #Correlation across genes DE in either spp.

dr.ph$Sig<-factor(paste(dr.ph$sig.DeRo,dr.ph$sig.PhHa))#create factor to colour DE genes by spp.
levels(dr.ph$Sig)<-c("None","PhHa","DeRo","Both")
g.dr.ph<-ggplot(dr.ph,aes(x=log2FoldChange.DeRo,y=log2FoldChange.PhHa))+
  cust.theme()+
  geom_abline(intercept=0,slope=1,linewidth=0.5,linetype='dotted')+
  geom_hline(yintercept=0,linewidth=0.5)+geom_vline(xintercept=0,linewidth=0.25)+
  geom_point(size=1,aes(colour=Sig))+
  scale_colour_manual(values=c('#aaaaaa',"#64e3dd","#e38888","#b383de"))+
  labs(caption=paste0("Comparison of ",nrow(dr.ph)," shared gene annotations"))+
  geom_text_repel(data=dr.ph[dr.ph$sig.DeRo & dr.ph$sig.PhHa,],
                  aes(x=log2FoldChange.DeRo,y=log2FoldChange.PhHa,label=gene),size=3,colour="#b383de")+
  xlim(c(-5,10))+ylim(c(-5,10))+annotate("text",label="Pearson's r = 0.53",x=7.5,y=-4)

ggsave('DeRo_PhHa_LPS_comp.png',dpi=600,height=12,width=10,
       plot=(g.dr.ph)/
         (plot_volcano(res.th.vs.c0)+ggtitle("DeRo")+plot_volcano(ph.th.vs.c0)+ggtitle("PhHa"))/
         (plot_GO_overrep_upreg(res.th.vs.c0)+plot_GO_overrep_upreg(ph.th.vs.c0))+
         plot_layout(heights=c(1.25,0.6,0.6)))


#AGE - C0
dds <- DESeqDataSetFromMatrix(
  countData = round(cts)[,samples$Trt=="C0" & samples$Sex=="F"],
  colData = samples[samples$Trt=="C0" & samples$Sex=="F",],
  design= ~ Est.Age)
dds <- DESeq(dds)
resultsNames(dds)
summary((results(dds, name="Est.Age",alpha=0.05)))
res.c0.age<-(data.frame(results(dds, name="Est.Age",alpha=0.05)))
res.c0.age$gene<-rownames(res.c0.age)
res.c0.age$sig<-res.c0.age$padj<0.05

write.csv(res.c0.age,"DeRo_C0_age.csv",row.names=FALSE,quote=FALSE)

#AGE - LPS
dds <- DESeqDataSetFromMatrix(
  countData = round(cts)[,samples$Trt=="LPS" & samples$Sex=="F"],
  colData = samples[samples$Trt=="LPS" & samples$Sex=="F",],
  design= ~ Est.Age)
dds <- DESeq(dds)
resultsNames(dds)
summary((results(dds, name="Est.Age",alpha=0.05)))
res.LPS.age<-(data.frame(results(dds, name="Est.Age",alpha=0.05)))
res.LPS.age$gene<-rownames(res.LPS.age)
res.LPS.age$sig<-res.LPS.age$padj<0.05
write.csv(res.LPS.age,"DeRo_LPS_age.csv",row.names=FALSE,quote=FALSE)

plot_volcano(res.c0.age)+ggtitle("C0 age")+plot_volcano(res.LPS.age)+ggtitle("LPS age")
ggsave("DeRo_age_volcanos.png",dpi=600,height=5,width=10)

samples.f<-samples[samples$Sex=="F",]
samples.f$dab1<-plotCounts(dds,gene="DAB1",intgroup="Trt",returnData = TRUE)$count
ggplot(samples.f,aes(x=Est.Age,y=dab1))+cust.theme()+geom_point()

#age-associated gene comp
ph.LPS.age<-read.csv("../nih_rna/july_25/results_files/phases3-7_TH_age_effect.csv",h=T)
dr.ph.LPS.age<-merge(res.LPS.age,ph.LPS.age,by='gene')
colnames(dr.ph.LPS.age)<-gsub("\\.x",".DeRo",colnames(dr.ph.LPS.age))
colnames(dr.ph.LPS.age)<-gsub("\\.y",".PhHa",colnames(dr.ph.LPS.age))
ggplot(dr.ph.LPS.age,aes(x=log2FoldChange.DeRo,y=log2FoldChange.PhHa))+cust.theme()+
  geom_point(aes(colour=paste(sig.DeRo,sig.PhHa)))+
  geom_text_repel(aes(label=gene))

cor.test(dr.ph.LPS.age$log2FoldChange.DeRo,dr.ph.LPS.age$log2FoldChange.PhHa,method="spearman")
