# Differential expression (DE) analysis of RNA-seq data from untreated and LPS-treated blood of Phyllostomus hastatus

# load_data_functions_summarise.R
## 00. load packages, define functions
## 01. read, filter, organise data
## 02. create summary plots

# run_DE_analyses.R
## 03. run DE analysis of paired untreated/LPS-treated samples
## 04. Spp comparison of LPS effect
## 05. Sex-specific LPS effect
## 06. Sex and age patterns in untreated samples
## 07. Sex and age patterns in LPS treated samples
## 08. Sex-specific age patterns in LPS treated samples

# run_WGCNA.R
## 09. WGCNA module analyses of LPS-treated samples


#####################

setwd("~/Documents/UMD_new/nih_rna/writeup/scripts_results_files/")

## load libraries and set working directory
library(DESeq2) 
library(edgeR) 
library(ggplot2) 
library(reshape2) 
library(patchwork)
library(gridExtra)
library(dplyr)
library(variancePartition)
library(sva)
library(clusterProfiler)
organism<-'org.Hs.eg.db'
library(organism, character.only = TRUE)
library(ggrepel)
library(ggridges)
library(ggbeeswarm)
library(ggplotify)
library(pheatmap)
library(ggupset)
library(tidyverse, warn.conflicts = FALSE)

## create functions and global variables for use in analysis
alpha=0.05

### add extra info to results files
gene_desc <- read.delim("gene_descriptions.tsv", header = TRUE, fill = TRUE)
gene_desc <- gene_desc[!duplicated(gene_desc$gene),]
colnames(gene_desc)[2]<-"description_Phha"
gene_desc$gene_symbol_human<-NA
gene_desc$description_human<-NA
recip.blastp<-read.csv("Hs_Ph_LOCs_reciprocalbesthits.csv",h=T)
colnames(recip.blastp)[3]<-'gene'
recip.blastp<-recip.blastp[,-c(1,2,7,8)]
gene_desc<-gene_desc[!gene_desc$gene %in% recip.blastp$gene,]
gene_desc<-rbind(recip.blastp,gene_desc)
add_gene_info<-function(x){
  #x<-x[!is.na(x$padj),]#get rid of NAs
  x$sig<-FALSE
  x[x$padj<alpha & !is.na(x$padj),]$sig<-TRUE
  x$gene<-rownames(x)
  return(merge(x,gene_desc,by='gene'))
}

### reduce legend sizes for GO overrepreesntation plots
addSmallLegend <- function(myPlot, pointSize = 1, textSize = 9, spaceLegend = 0.75) {
  myPlot + guides(shape = guide_legend(override.aes = list(size = pointSize)),
                  color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

### Define theme for plots
SexPalette<-c("#C74955","#5090bf")#Sex
TrtPalette<-c("#81ad50","#c987e8")#Treatment
StPalette<-c("#7eb2d9","#195f91")#Male status
cust.theme<-function(){
  theme_minimal() %+replace%
    theme(panel.grid=element_blank(),axis.ticks=element_line(),
          plot.background=element_rect(fill='white',colour='white'),
          panel.background=element_rect(fill="#fcfbfa",colour='black',linewidth = 0.25))
}

### main GO overrepresentation plot
charlim=45 #limit num. characters to show on X-axis (GO term descriptions)
plot_GO_fun<-function(go.OR){
  go.OR<-go.OR[order(go.OR$p.adjust),]
  go.OR<-go.OR[!duplicated(go.OR$geneID),]
  go.OR<-go.OR[order(go.OR$p.adjust),] %>% head(n=10)
  min1<-min(-log10(go.OR$p.adjust))*0.95
  max1<-max(-log10(go.OR$p.adjust))
  min2<-min((go.OR$Count))
  max2<-max((go.OR$Count))
  if(any(nchar(go.OR$Description)>charlim)){
    go.OR[nchar(go.OR$Description)>charlim,]$Description<-
      paste0(substr(go.OR[nchar(go.OR$Description)>charlim,]$Description,start=0,stop = charlim),"...")
  }
  return(
    addSmallLegend(
      ggplot(go.OR,aes(x=-log10(p.adjust),y=Description,fill=-log10(p.adjust)))+
        theme_minimal()+theme(axis.title.y=element_blank(),panel.grid=element_blank(),plot.background=element_rect(fill='white',colour='white'),
                              panel.background=element_rect(fill="white",colour='white',linewidth = 0.25),panel.border = element_blank())+
        geom_segment(aes(x=min1,xend=-log10(p.adjust),y=Description,group=Description),colour="#999",linewidth=0.5)+
        geom_point(aes(size=Count),shape=21)+
        scale_fill_gradient(low = "#eddd8e", high = "#ba3636", na.value = NA)+
        #scale_fill_viridis(option="A")+
        scale_size(range=c(3,7),limits=c(min2,max2))+
        # theme(axis.title.y=element_blank(),legend.title=element_text(size=10),axis.ticks.y=element_line(),
        #       legend.position='right',panel.grid=element_blank(),plot.background=element_rect(fill='white',colour='white'),
        #       panel.background=element_rect(fill="#fcfbfa",colour='black',linewidth = 0.25),axis.text.y=element_text(size=8.5),
        #       legend.frame = element_rect(color = "#555555", size = 0.2))+
        xlab("-log10(Padj)")+ guides(fill = "none")
    )
  )
}


### run GO analysis 
plot_GO_overrep<-function(results_table,dir){
  gene_list_interest<-results_table[results_table$padj<alpha & !is.na(results_table$padj),]
  if(dir=="up"){
    gene_list_interest<-gene_list_interest[gene_list_interest$log2FoldChange>0,]$gene
  }else if (dir=="dn"){
    gene_list_interest<-gene_list_interest[gene_list_interest$log2FoldChange<0,]$gene
  }
  go.OR<-clusterProfiler::simplify(enrichGO(gene = gene_list_interest,
                                            universe = results_table$gene,#list of all genes
                                            keyType = "SYMBOL",
                                            OrgDb = organism,
                                            ont = "BP",
                                            pAdjustMethod = "BH",
                                            pvalueCutoff = 1,
                                            readable = TRUE),
                                   cutoff=0.7,by = "p.adjust",select_fun = min,
                                   measure = "Wang",semData = NULL)@result
  plot_GO_fun(go.OR)
}

### plot volcano
plot_volcano<-function(df){
  if (any(df$pvalue==0 & !is.na(df$pvalue))){
    df[df$pvalue==0 & !is.na(df$pvalue),]$pvalue<-1.520048e-299#replace 0 values for finite log-transformation
  }
  df<-df[!is.na(df$padj),]
  anno.set<-df[abs(df$log2FoldChange)>quantile(abs(df$log2FoldChange),0.95) & df$padj<0.05 | 
                 -log10(df$pvalue)>quantile(-log10(df$pvalue),0.95)  & df$padj<0.05,]
  anno.set<-anno.set[-grep("^LOC",anno.set$gene),]#don't bother display non-functionally annotated genes
  g.volcano<-ggplot(df,aes(x=log2FoldChange,y=-log10(pvalue)))+
    cust.theme()+theme(legend.position='none')+
    geom_vline(xintercept=0,linetype='dashed')+
    geom_point(aes(colour=padj<alpha),size=1.5,alpha=0.75)+scale_colour_manual(values=c("#aaaaaa","#e07575"))+
    geom_text_repel(data=anno.set,colour='black',fontface = 'italic',
                    aes(label=gene),size=2.5,max.overlaps = 8,min.segment.length = 0)
  return(g.volcano)
}

# 01. read, filter, organise data
#####################
# 01. READ. FILTER, ORGANISE DATA

samples.all<-read.csv("samples_all.csv",h=T,row.names=1)
cts.all<-read.csv("cts_all.csv",h=T,row.names=1)

## scale age, ensure variables are treated appropriately (e.g., band and year treated as factors)
samples.all$Age<-scale(samples.all$Est.Age)
samples.all$Band<-(factor(samples.all$Band))
samples.all$Yr<-(factor(samples.all$Yr))
colnames(samples.all)[which(colnames(samples.all)=="LN.N.L..")]<-"lnNLR"
samples.all$Phase<-samples.all$Phase.batch
samples.all$Batch<-samples.all$Phase
samples.all$ID<-samples.all$Band

## make subsets
cts.p7<-cts.all[,samples.all$Phase=="Phase7.1" | samples.all$Phase=="Phase7.2"]
samples.p7<-samples.all[samples.all$Phase=="Phase7.1"| samples.all$Phase=="Phase7.2",]
cts.all.th<-cts.all[,samples.all$Trtmt=="TH"]
samples.all.th<-samples.all[samples.all$Trtmt=="TH",]
cts.all.c0<-cts.all[,samples.all$Trtmt=="C0"]
samples.all.c0<-samples.all[samples.all$Trtmt=="C0",]
cts.all<-cts.all[,samples.all$Trtmt %in% c("C0","TH")]
samples.all<-samples.all[samples.all$Trtmt %in% c("C0","TH"),]

#####################
# 02. SUMMARY PLOTS

## plot age distributions
g.ages.th<-ggplot(samples.all.th,aes(x=Est.Age,fill=Sex))+cust.theme()+geom_histogram(alpha=0.75,colour='black')+
  scale_fill_manual(values=SexPalette)+facet_grid(Sex~.)+theme(legend.position='none')+
  ggtitle("LPS-treated samples")
g.ages.c0<-ggplot(samples.all.c0,aes(x=Est.Age,fill=Sex))+cust.theme()+geom_histogram(alpha=0.75,colour='black')+
  scale_fill_manual(values=SexPalette)+facet_grid(Sex~.)+theme(legend.position='none')+
  ggtitle("Untreated samples")
g.ages.paired<-ggplot(samples.all[samples.all$Phase %in% c("Phase7.1","Phase7.2") & samples.all$Trtmt=="C0",],aes(x=Est.Age,fill=Sex))+
  cust.theme()+geom_histogram(alpha=0.75,colour='black')+
  scale_fill_manual(values=SexPalette)+facet_grid(Sex~.)+theme(legend.position='none')+
  ggtitle("Paired individuals")

g.ages.th+g.ages.c0+g.ages.paired+plot_layout(nrow=3)
#ggsave("~/Documents/UMD_new/nih_rna/writeup/sex_age_distribution.png",width=4,height=8)

## variance partition - exploratory (takes a while to run)

# form <- ~  Age + lnNLR + WT + (1 |Phase.batch) + (1 | Sex) + (1 | Cave)
# varPart <- fitExtractVarPartModel(cts.all.th, form, samples.all.th,n=1000)
# vp <- sortCols(varPart)
# vp<-plotVarPart(vp)$data
# g2<-ggplot(vp,aes(x=variable,y=value,fill=variable))+
#   theme_bw()+theme(panel.grid=element_blank())+
#   theme(axis.text.x=element_text(angle = 45,vjust=1, hjust=1,size=8),axis.title.y=element_blank(),
#         panel.background = element_rect(fill='white',colour='#aaaaaa',linewidth=0.25))+
#   geom_quasirandom(size=0.5,alpha=0.5,colour='#888888')+
#   geom_boxplot(colour='black',width=0.5,alpha=0.75,outlier.shape=NA)+
#   theme(legend.position='none')+ylab("Variance explained (%)")
# g2<-g2+coord_flip(ylim=c(0,100)) 
