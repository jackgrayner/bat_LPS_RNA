library(DESeq2)
library(dplyr)

add_gene_info<-function(x){
  x<-data.frame(x)
  x$sig<-FALSE
  if(nrow(x[x$padj<0.05 & !is.na(x$padj),]>0)){
    x[x$padj<0.05 & !is.na(x$padj),]$sig<-TRUE
  }
  x$miRNA<-gsub("_.*","",rownames(x))
  x$spp<-substr(x$miRNA,0,3)
  x<-x[order(x$pvalue),]
  return(x)
}
cust.theme<-function(){
  theme_minimal() %+replace%
    theme(legend.position='none',plot.title=element_text(size=14,hjust = 0,vjust = 0.5),
          plot.background=element_rect(fill='white',colour='white'),
          panel.background=element_rect(fill="#fcfbfa",colour='black',linewidth = 0.25))
}

samples<-read.csv("./final_counts/sample_info.csv",h=T)
rownames(samples)<-samples$sample

cts<-read.table("./final_counts/all_samples_expression.tsv",h=T,sep="\t")

#sum cts across duplicate mature miRNAs
cts <- cts %>%
  group_by(miRNA) %>%
  summarise(across(starts_with("read_count"), sum))
cts<-data.frame(cts)

cts<-data.frame(cts,row.names=1)
colnames(cts)<-gsub("X","",colnames(cts))
cts<-cts[rowSums(cts)>20,]

summary(colnames(cts) == rownames(samples))

#DE analysis
dds <- DESeq(DESeqDataSetFromMatrix(countData = cts,
                                    colData = samples,
                                    design= ~ trt))
resultsNames(dds)
summary(results(dds, name="trt_TH_vs_C0",alpha=0.05))
res.trt<- add_gene_info(results(dds, name="trt_TH_vs_C0",alpha=0.05))

pca.trt<-plotPCA(vst(dds,nsub=200),intgroup="trt",returnData=TRUE)
pca.trt$sex<-samples$Sex
pca.trt$age<-samples$Est..Age
pca.trt$trt<-factor(pca.trt$trt)
levels(pca.trt$trt)=c("none","LPS")
g.pca<-ggplot(pca.trt,aes(x=PC1,y=PC2,colour=trt))+
  cust.theme()+geom_point(size=2)+scale_colour_manual(values=TrtPalette)+
  xlab("PC1 (19% var.)")+ylab("PC2 (15% var,)")+
  theme(legend.position='right')+labs(colour="Treatment")

anno.set<-res.trt[res.trt$padj<0.05,]
g.volcano<-ggplot(res.trt[!is.na(res.trt$padj),],aes(x=log2FoldChange,y=-log10(pvalue)))+cust.theme()+
  geom_vline(xintercept=0,linetype='dashed',colour="#aaaaaa")+
  geom_point(aes(colour=padj<0.05),size=1.5,alpha=1)+scale_colour_manual(values=c("#aaaaaa","#C74955"))+
  geom_text_repel(data=anno.set,colour='black',fontface = 'italic',
                  aes(label=miRNA),size=2,max.overlaps = 20,min.segment.length = 0,
                  segment.colour = "#aaaaaa",segment.size = 0.5,force_pull = 0.5)

g.pca+labs(tag="A")+g.volcano+labs(tag="B")
ggsave('volcano_PCA.png',dpi=600,height=3,width=6.5)



View(res.trt)
write.csv(res.trt,"./results_files/miRNA_trt_effect.csv",quote=FALSE)

# test sex/age effect
dds <- DESeq(DESeqDataSetFromMatrix(countData = cts[,samples$trt=="TH"],
                                    colData = samples[samples$trt=="TH",],
                                    design= ~ Sex + Est..Age #+ lnNLR
                                      ))
resultsNames(dds)
summary(results(dds, name="Sex_M_vs_F",alpha=0.05))
summary(results(dds, name="Est..Age",alpha=0.05))
summary(results(dds, name="lnNLR",alpha=0.05))
res.th.age<- add_gene_info(results(dds, name="Est..Age",alpha=0.05))
res.th.age[res.th.age$padj<0.05 & !is.na(res.th.age$padj),]
res.th.age[res.th.age$sig,]
write.csv(res.th.age,"./results_files/miRNA_TH_samples_age_effect.csv",quote=FALSE)

g.volcano<-ggplot(res.th.age[!is.na(res.th.age$padj),],aes(x=log2FoldChange,y=-log10(pvalue)))+cust.theme()+
  geom_vline(xintercept=0,linetype='dashed')+
  geom_point(aes(colour=padj<0.05),size=1.5,alpha=1)+scale_colour_manual(values=c("#aaaaaa","#C74955"))+
  geom_text_repel(colour='black',fontface = 'italic',
                  aes(label=miRNA),size=2,max.overlaps = 10,min.segment.length = 0,
                  segment.colour = "#aaaaaa",segment.size = 0.5,force_pull = 0.5)

vst<-assay(vst(dds,nsub=200,blind=FALSE))
vst.mir7<-data.frame(count=(vst[rownames(vst)=="hsa-miR-7-5p",]),sex=samples[samples$trt=="TH",]$Sex,age=samples[samples$trt=="TH",]$Est..Age)
vst.mir599<-data.frame(count=(vst[rownames(vst)=="pal-miR-599-3p",]),sex=samples[samples$trt=="TH",]$Sex,age=samples[samples$trt=="TH",]$Est..Age)
g.mir7<-ggplot(vst.mir7,aes(x=age,y=count,colour=sex))+cust.theme()+
  geom_point()+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+
  ggtitle("hsa-miR-7-5p")
g.mir599<-ggplot(vst.mir599,aes(x=age,y=count,colour=sex))+cust.theme()+
  geom_point()+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+
  ggtitle("pal-miR-599-3p")

g.mir7+g.mir599
ggsave('age_assoc_mirnas.png',dpi=600,height=3,width=6)

