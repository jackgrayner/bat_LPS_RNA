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

samples<-read.csv("~/Documents/UMD_new/nih_rna/miRNA/sample_info.csv",h=T)
rownames(samples)<-samples$sample

cts<-read.table("~/Documents/UMD_new/nih_rna/miRNA/all_samples_expression.tsv",h=T,sep="\t")

#sum cts across duplicate mature miRNAs
cts <- cts %>%
  group_by(miRNA) %>%
  summarise(across(starts_with("read_count"), sum))
cts<-data.frame(cts)

rownames(cts)<-cts$miRNA
cts<-cts[,-1]
colnames(cts)<-rownames(samples)
cts<-cts[rowSums(cts)>20,]
cts<-cts[order(rowSums(cts),decreasing = TRUE),]

head(cts)

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


res.trt[res.trt$sig,]
summary(factor(res.trt[res.trt$sig,]$spp))
summary(factor(res.trt$spp))
res.trt$bat_annotated<-res.trt$spp %in% c("pal","efu","aja","mmy","mlu")
table(res.trt$sig,res.trt$bat)
fisher.test(table(res.trt$sig,res.trt$bat_annotated))
View(res.trt[grep("hsa",res.trt$miRNA),])
hsa.sig<-res.trt[grep("hsa",res.trt$miRNA),] %>% filter(sig)
write.table(hsa.sig$miRNA,"~/Documents/UMD_new/nih_rna/miRNA/hsa_sig.txt",quote=FALSE,row.names = FALSE,col.names = FALSE)

anno.set<-res.trt[res.trt$padj<0.05,]
g.volcano<-ggplot(res.trt[!is.na(res.trt$padj),],aes(x=log2FoldChange,y=-log10(pvalue)))+cust.theme()+
  geom_vline(xintercept=0,linetype='dashed',colour="#aaaaaa")+
  geom_point(aes(colour=padj<0.05),size=1.5,alpha=1)+scale_colour_manual(values=c("#aaaaaa","#C74955"))+
  geom_text_repel(data=anno.set,colour='black',fontface = 'italic',
                  aes(label=miRNA),size=2,max.overlaps = 20,min.segment.length = 0,
                  segment.colour = "#aaaaaa",segment.size = 0.5,force_pull = 0.5)

g.pca+labs(tag="A")+g.volcano+labs(tag="B")
ggsave('~/Documents/UMD_new/nih_rna/miRNA/volcano_PCA.png',dpi=600,height=3,width=6.5)



View(res.trt)
write.csv(res.trt,"~/Documents/UMD_new/nih_rna/miRNA/miRNA_trt_effect.csv",quote=FALSE)


#human pro/anti-inflammatory mirnas
pros<-c("hsa-miR-155","hsa-miR-92a","hsa-miR-200","hsa-miR-23a",
        "hsa-miR-27a","hsa-miR-29c","hsa-miR-138","hsa-miR-34a",
        "hsa-miR-34c","hsa-miR-132","hsa-let-7a",
        "mmu-miR-155","mmu-miR-92a","mmu-miR-200","mmu-miR-23a",
        "mmu-miR-27a","mmu-miR-29c","mmu-miR-138","mmu-miR-34a",
        "mmu-miR-34c","mmu-miR-132","mmu-let-7a",
        "rno-miR-155","rno-miR-92a","rno-miR-200","rno-miR-23a",
        "rno-miR-27a","rno-miR-29c","rno-miR-138","rno-miR-34a",
        "rno-miR-34c","rno-miR-132","rno-let-7a")

antis<-c("hsa-miR-10a","hsa-miR-7","hsa-miR-126","hsa-miR-146a",
        "hsa-miR-124","hsa-miR-125b","hsa-miR-31","hsa-miR-210",
        "hsa-miR-24","hsa-miR-149","hsa-miR-181","hsa-miR-150",
        "hsa-miR-143","hsa-miR-9","hsa-miR-142","hsa-miR-223","hsa-miR-21",
        "mmu-miR-10a","mmu-miR-7","mmu-miR-126","mmu-miR-146a",
        "mmu-miR-124","mmu-miR-125b","mmu-miR-31","mmu-miR-210",
        "mmu-miR-24","mmu-miR-149","mmu-miR-181","mmu-miR-150",
        "mmu-miR-143","mmu-miR-9","mmu-miR-142","mmu-miR-223","mmu-miR-21",
        "rno-miR-10a","rno-miR-7","rno-miR-126","rno-miR-146a",
        "rno-miR-124","rno-miR-125b","rno-miR-31","rno-miR-210",
        "rno-miR-24","rno-miR-149","rno-miR-181","rno-miR-150",
        "rno-miR-143","rno-miR-9","rno-miR-142","rno-miR-223","rno-miR-21")

res.trt$inflam<-NA
res.trt[gsub("-[0-9]p","",res.trt$miRNA) %in% pros,]$inflam<-"pro_inflam"
res.trt[gsub("-[0-9]p","",res.trt$miRNA) %in% antis,]$inflam<-"anti_inflam"
g.volcano<-ggplot(res.trt[!is.na(res.trt$inflam),],aes(x=log2FoldChange,y=-log10(padj)))+cust.theme()+
  geom_vline(xintercept=0,linetype='dashed')+theme(legend.position = 'right')+
  geom_hline(yintercept=1.3,linetype='dashed')+
  #geom_point(aes(colour=padj<0.05),size=1,alpha=1)+scale_colour_manual(values=c("#aaaaaa","darkred"))+
  geom_point(aes(x=log2FoldChange,y=-log10(padj),colour=inflam),size=2,alpha=1)+
  geom_text_repel(colour='black',fontface = 'italic',
                  aes(label=miRNA),size=3,max.overlaps = 1,min.segment.length = 0)


#human macrophage comparison
hs.1hr<-read.csv("~/Documents/UMD_new/nih_rna/rodents/miRNA_hs_macro_1hr.csv")
hs.8hr<-read.csv("~/Documents/UMD_new/nih_rna/rodents/miRNA_hs_macro_8hr.csv")
head(hs.1hr)
head(res.trt)
summary(res.trt$miRNA %in% hs.1hr$miRNA)

#

dds <- DESeq(DESeqDataSetFromMatrix(countData = cts[,samples$trt=="TH"],
                                    colData = samples[samples$trt=="TH",],
                                    design= ~ Sex + Est..Age + lnNLR
                                      ))
resultsNames(dds)
summary(results(dds, name="Sex_M_vs_F",alpha=0.05))
summary(results(dds, name="Est..Age",alpha=0.05))
summary(results(dds, name="lnNLR",alpha=0.05))
res.th.age<- add_gene_info(results(dds, name="Est..Age",alpha=0.05))
res.th.age[res.th.age$padj<0.05 & !is.na(res.th.age$padj),]
res.th.age[res.th.age$sig,]
write.csv(res.th.age,"~/Documents/UMD_new/nih_rna/miRNA/miRNA_TH_samples_age_effect.csv",quote=FALSE)

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
ggsave('~/Documents/UMD_new/nih_rna/miRNA/age_assoc_mirnas.png',dpi=600,height=3,width=6)

lm.vst<-lm(count~sex+age,data=vst)
car::Anova(lm.vst,type="II")
