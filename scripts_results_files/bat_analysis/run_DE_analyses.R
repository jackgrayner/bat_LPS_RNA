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

#####################
# 03. DE ANALYSIS OF PAIRED SAMPLES

## Run model
dds.p7 <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.p7),
  colData = samples.p7,
  design= ~ Trtmt+Band))
resultsNames(dds.p7)
summary(results(dds.p7, name="Trtmt_TH_vs_C0",alpha=alpha))
res.th.vs.c0<- add_gene_info(data.frame(results(dds.p7, name="Trtmt_TH_vs_C0",alpha=alpha)))
res.th.vs.c0<-res.th.vs.c0[order(res.th.vs.c0$pvalue),]
lps.genes<-res.th.vs.c0[res.th.vs.c0$sig,]$gene

## save results and vst counts
write.csv(res.th.vs.c0,"phase7_paired_th_vs_c0.csv",row.names = FALSE)
write.csv(assay(vst(dds.p7,blind=FALSE)),"phase7_paired_vst.csv",row.names = FALSE)

## inspecting some genes of interest
res.th.vs.c0[res.th.vs.c0$gene %in% c("NFKB1","IL1A","IL1B","IL6"),]
res.th.vs.c0[res.th.vs.c0$gene %in% c("IL23R","TNFRSF11A","SLC11A1"),]

## plot PCA
pca.trt<-plotPCA(vst(dds.p7),intgroup="Trtmt",returnData=TRUE,ntop=1000)
pca.trt$sex<-samples.p7$Sex
pca.trt$age<-samples.p7$Est.Age
pca.trt$Band<-samples.p7$Band
pca.trt$Phase.batch<-samples.p7$Phase.batch
pca.trt<-pca.trt[pca.trt$Trtmt %in% c("C0","TH"),]
pca.trt$trt<-factor(pca.trt$Trtmt)
levels(pca.trt$trt)=c("none","LPS")
p7.pca<-ggplot(pca.trt,aes(x=PC1,y=PC2,fill=trt))+
  cust.theme()+
  geom_line(aes(group=Band),colour="#aaaaaa",linewidth = 0.25)+
  geom_point(size=2.75,shape=21)+scale_fill_manual(values=TrtPalette)+
  xlab("PC1 (51% var.)")+ylab("PC2 (8% var.)")+
  theme(legend.position='left')+labs(colour="Treatment")

## plot differences in PC1 between untreated and LPS treated samples
pca.paired<-merge(pca.trt[pca.trt$Trtmt=="C0",],pca.trt[pca.trt$Trtmt=="TH",],by="Band")
pca.paired$diff<-pca.paired$PC1.y-pca.paired$PC1.x
p7.pcdiff<-ggplot(pca.paired,aes(x=age.x,y=diff,colour=sex.x))+
  cust.theme()+geom_point()+geom_smooth(method='lm')+
  scale_colour_manual(values=SexPalette)+
  labs(colour="Sex",x="Age",y="PC1 difference (LPS - untreated)")

p7.pca+labs(tag="A")+p7.pcdiff+labs(tag="B")
#ggsave("~/Documents/UMD_new/nih_rna/writeup/paired_PC1_diff_sexage.png",dpi=600,height=3.5,width=8)

## run LM of PC1, testing sex and age effects
lm.pca.paired<-lmer(diff~sex.x*age.x+(1|Phase.batch.x),data=pca.paired)
plot(lm.pca.paired)
car::Anova(lm.pca.paired,type="III")#borderline sig. sex * age interaction

### Plot Volcano  for TH-affected genes
g.volcano<-plot_volcano(res.th.vs.c0)

### Plot GO overrepresentation for TH-affected genes
go.all.th.up<-plot_GO_overrep(res.th.vs.c0,'up')
#plot_GO_overrep(res.th.vs.c0,'dn')

### run and plot KEGG enrichment for TH-affected genes
ids<-bitr(res.th.vs.c0[res.th.vs.c0$padj<alpha & res.th.vs.c0$log2FoldChange>0,]$gene,
          fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism) %>% filter(!duplicated(ENTREZID))
df2 = res.th.vs.c0[res.th.vs.c0$gene %in% ids$SYMBOL,] %>% mutate(SYMBOL=gene)
df2<-merge(df2,ids,by='SYMBOL')
df2<-df2[order(df2$log2FoldChange,decreasing = TRUE),]
kegg_gene_list <- df2$log2FoldChange
names(kegg_gene_list) <- df2$ENTREZID
geneList = sort(kegg_gene_list, decreasing = TRUE)
kk2 <- gseKEGG(geneList     = (geneList),
               organism     = 'human',
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid")@result
kk2<-head(kk2,n=10)
kk2$Count<-kk2$setSize
### set axis lims 
min1<-min(-log10(kk2$p.adjust))-0.5
max1<-max(-log10(kk2$p.adjust))+0.5
min2<-min((kk2$setSize))
max2<-max((kk2$setSize))
dotplot.kegg<-ggplot(kk2,aes(x=-log10(p.adjust),y=Description,fill=-log10(p.adjust)))+
  theme_minimal()+theme(axis.title.y=element_blank(),panel.grid=element_blank(),plot.background=element_rect(fill='white',colour='white'),
                        panel.background=element_rect(fill="white",colour='white',linewidth = 0.25),panel.border = element_blank())+
  geom_segment(aes(x=min1,xend=-log10(p.adjust),y=Description,group=Description),colour="#999",linewidth=0.5)+
  geom_point(aes(size=Count),shape=21)+
  scale_fill_gradient(low = "#eddd8e", high = "#ba3636", na.value = NA)+
  scale_size(range=c(3,7),limits=c(min2,max2))+
  xlim(c(min1,max1))+xlab("-log10(Padj)")+ guides(fill = "none")

#####################
# 03.5. PCA of all samples
## Run model
dds.all <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all),
  colData = samples.all,
  design= ~ 1))
resultsNames(dds.all)

pca.trt<-plotPCA(vst(dds.all),intgroup="Trtmt",returnData=TRUE,ntop=1000)
#plotPCA(vst(dds.all),intgroup="Trtmt",returnData=FALSE,ntop=1000) #to get %var explained
pca.trt$sex<-samples.all$Sex
pca.trt$age<-samples.all$Est.Age
pca.trt$Band<-samples.all$Band
pca.trt$Phase.batch<-samples.all$Phase.batch
pca.trt$lnNLR<-samples.all$lnNLR
pca.trt$status<-samples.all$STATUS
pca.trt<-pca.trt[pca.trt$Trtmt %in% c("C0","TH"),]
pca.trt$trt<-factor(pca.trt$Trtmt)
levels(pca.trt$trt)=c("none","LPS")

all.pca<-ggplot(pca.trt,aes(x=PC1,y=PC2,fill=trt))+
  cust.theme()+
  geom_line(aes(group=Band),colour="#aaaaaa",linewidth = 0.25)+
  geom_point(size=2.75,shape=21)+scale_fill_manual(values=TrtPalette)+
  xlab("PC1 (47% var.)")+ylab("PC2 (8% var.)")+
  theme(legend.position='left')+labs(colour="Treatment")

#lm.pca.all<-lmer(PC1~sex*age+(1|Phase.batch),data=pca.trt[pca.trt$Trtmt=="TH",])
lm.pca.all<-lmer(PC1~sex+age+(1|Phase.batch),data=pca.trt[pca.trt$Trtmt=="TH",])
lm.pca.all<-lmer(PC1~sex+age+lnNLR+(1|Phase.batch),data=pca.trt[pca.trt$Trtmt=="TH",])
plot(lm.pca.paired)
qqnorm(resid(lm.pca.paired))
qqline(resid(lm.pca.paired))
#car::Anova(lm.pca.all,type="III")#no sex * age interaction
car::Anova(lm.pca.all,type="II")


#####################
# 04. SPP COMPARISON -> create multi spp heatmap of LPS genes

lps<-read.table("lps_GO_genes.txt",h=F)
lps<-lps$V1

## bat LPS results come from above model, others from independent analyses of existing data
bat.lfc<-res.th.vs.c0
bab.lfc<-read.csv("./baboon_reanalysis/baboon_LPS.csv") %>% mutate(log2FoldChange=-log2FoldChange)#correct contrast direction
mac.lfc<-read.csv("./macaque_reanalysis/macaque_LPS.csv") %>% mutate(log2FoldChange=-log2FoldChange)
pig.lfc<-read.csv("./pigs_reanalysis/pig_LPS.csv")

## define significant genes in each spp.
bat.sig<-res.th.vs.c0[res.th.vs.c0$padj<0.05  & !is.na(res.th.vs.c0$padj),]$gene
bab.sig<-bab.lfc[bab.lfc$padj<0.05 & !is.na(bab.lfc$padj),]$gene
mac.sig<-mac.lfc[mac.lfc$padj<0.05 & !is.na(mac.lfc$padj),]$gene
pig.sig<-pig.lfc[pig.lfc$padj<0.05 & !is.na(pig.lfc$padj),]$gene

## create vector of LPS GO term genes sig. in at least one spp. (or in multiple?)
lps<-lps[lps %in% c(bat.sig,bab.sig,mac.sig,pig.sig)]

## create a heatmap of the LPS GO term genes
bat.lfc1<-bat.lfc[bat.lfc$gene %in% lps,c("gene","log2FoldChange")] %>% `colnames<-`(c("gene","bat"))
bab.lfc1<-bab.lfc[bab.lfc$gene %in% lps,c("gene","log2FoldChange")] %>% `colnames<-`(c("gene","baboon"))
mac.lfc1<-mac.lfc[mac.lfc$gene %in% lps,c("gene","log2FoldChange")] %>% `colnames<-`(c("gene","macaque"))
pig.lfc1<-pig.lfc[pig.lfc$gene %in% lps,c("gene","log2FoldChange")] %>% `colnames<-`(c("gene","pig"))

imp.gene.lfc<-left_join(bat.lfc1,pig.lfc1,by="gene") %>% left_join(.,bab.lfc1,by='gene') %>% left_join(.,mac.lfc1,by='gene')
imp.gene.lfc <- data.frame(imp.gene.lfc[,-1], row.names = imp.gene.lfc[,1])
imp.gene.lfc<-imp.gene.lfc[complete.cases(imp.gene.lfc),]

neg_colors <- colorRampPalette(c("#5a9fdb", "#ffffff"))(10)
pos_colors <- colorRampPalette(c("#f7f2e1","#f2d479","#f5ad1d", "#d62302"))(20)
my_colors  <- c(neg_colors, pos_colors)

breaks <- c(
  seq(min(imp.gene.lfc, na.rm = TRUE), 0, length.out = 10),
  seq(0, max(imp.gene.lfc, na.rm = TRUE), length.out = 20)[-1]
)

## transpose dataframe, add significant annotations
imp.gene.lfc.t<-data.frame(t(imp.gene.lfc))
anno_col.lps<-data.frame(lps.gene=factor(as.integer(colnames(imp.gene.lfc.t) %in% lps)))
rownames(anno_col.lps)<-colnames(imp.gene.lfc.t)
anno_col.dge<-data.frame(macaque=case_when(colnames(imp.gene.lfc.t) %in% mac.sig ~ "Y",.default = "N"),
                         baboon=case_when(colnames(imp.gene.lfc.t) %in% bab.sig ~ "Y",.default = "N"),
                         pig=case_when(colnames(imp.gene.lfc.t) %in% pig.sig ~ "Y",.default = "N"),
                         bat=case_when(colnames(imp.gene.lfc.t) %in% bat.sig ~ "Y",.default = "N")
)
rownames(anno_col.dge)<-colnames(imp.gene.lfc.t)
my_colour=list('macaque' = c("N"="white","Y"='black'),
               'baboon' = c("N"="white","Y"='black'),
               'pig' = c("N"="white","Y"='black'),
               'bat' = c("N"="white","Y"='black'))

spp.lps.heatmap<-(pheatmap(imp.gene.lfc.t,cluster_rows = FALSE,cluster_cols = TRUE,
                           border_color = "#666",treeheight_row = 10,angle_col = 90,treeheight_col = 0,cutree_cols = 9,
                           color = my_colors,breaks=breaks,na_col = "#bbbbbb",fontsize_row = 8,fontsize_col = 5.5,
                           annotation_col = anno_col.dge,annotation_colors = my_colour,annotation_legend = FALSE))

#ggsave('species_comparison_heatmap_hori.svg',height=1.5,width=9.25,plot=spp.heatmap.t)

## now create plots of correlations in LPS response across one-to-one orthologs

bat.lfc<-read.csv("./ortho_dge/bat_ortho_dge.csv")
bat.sig1<-bat.lfc[bat.lfc$sig,]$gene
bab.lfc<-read.csv("./ortho_dge/baboon_ortho_dge.csv")%>% mutate(log2FoldChange=-log2FoldChange)#correct contrast direction
bab.sig1<-bab.lfc[bab.lfc$sig,]$gene
mac.lfc<-read.csv("./ortho_dge/macaque_ortho_dge.csv") %>% mutate(log2FoldChange=-log2FoldChange)
mac.sig1<-mac.lfc[mac.lfc$sig,]$gene
pig.lfc<-read.csv("./ortho_dge/pig_ortho_dge.csv")
pig.sig1<-pig.lfc[pig.lfc$sig,]$gene

bat.lfc2<-bat.lfc[,c("gene","log2FoldChange","eggNOG_OGs")] %>% `colnames<-`(c("batgene","bat","eggNOG_OGs"))
bab.lfc2<-bab.lfc[,c("gene","log2FoldChange","eggNOG_OGs")] %>% `colnames<-`(c("babgene","baboon","eggNOG_OGs"))
mac.lfc2<-mac.lfc[,c("gene","log2FoldChange","eggNOG_OGs")] %>% `colnames<-`(c("macgene","macaque","eggNOG_OGs"))
pig.lfc2<-pig.lfc[,c("gene","log2FoldChange","eggNOG_OGs")] %>% `colnames<-`(c("piggene","pig","eggNOG_OGs"))

imp.gene.lfc<-left_join(bat.lfc2,pig.lfc2,by="eggNOG_OGs") %>% left_join(.,bab.lfc2,by='eggNOG_OGs') %>% left_join(.,mac.lfc2,by='eggNOG_OGs')
rownames(imp.gene.lfc)<-imp.gene.lfc$batgene

###  correlations - all
cor.test(imp.gene.lfc$bat,imp.gene.lfc$pig,method="spearman")
cor.test(imp.gene.lfc$bat,imp.gene.lfc$macaque,method="spearman")
cor.test(imp.gene.lfc$bat,imp.gene.lfc$baboon,method="spearman")

imp.gene.lfc$gene<-rownames(imp.gene.lfc)

### create pairwise data frames and plots
imp.gene.lfc$bat.pig.sig<-paste(imp.gene.lfc$batgene %in% bat.sig1,imp.gene.lfc$piggene %in% pig.sig1)
imp.gene.lfc$bat.macaque.sig<-paste(imp.gene.lfc$batgene %in% bat.sig1,imp.gene.lfc$macgene %in% mac.sig1)
imp.gene.lfc$bat.baboon.sig<-paste(imp.gene.lfc$batgene %in% bat.sig1,imp.gene.lfc$babgene %in% bab.sig1)

g.cor.pig<-ggplot(imp.gene.lfc,aes(x=bat,y=pig))+cust.theme()+xlab("bat LogFC")+ylab("pig LogFC")+
  geom_point(data=imp.gene.lfc[!imp.gene.lfc$bat.pig.sig=="FALSE FALSE",],aes(x=bat,y=pig,colour=bat.pig.sig),alpha=0.85,size=1.5)+
  geom_hline(yintercept=0,linetype='dotted')+geom_vline(xintercept=0,linetype='dotted')+
  theme(legend.position=c(4,2.5))+scale_colour_manual(values=c("#0C7C59","#444D5D","#EF6F6C"))#+#scale_color_manual(values=c("#888","#555"))+
#geom_text_repel(aes(label=gene),size=2,fontface='italic',alpha=1,max.overlaps = 10,segment.size = 0.125)

g.cor.mac<-ggplot(imp.gene.lfc,aes(x=bat,y=macaque))+cust.theme()+xlab("bat LogFC")+ylab("macaque LogFC")+
  geom_point(data=imp.gene.lfc[!imp.gene.lfc$bat.macaque.sig=="FALSE FALSE",],
             aes(x=bat,y=macaque,colour=bat.macaque.sig),alpha=0.85,size=1.5)+
  geom_hline(yintercept=0,linetype='dotted')+geom_vline(xintercept=0,linetype='dotted')+
  theme(legend.position=c(4,2.5))+scale_colour_manual(values=c("#0C7C59","#444D5D","#EF6F6C"))+
  ylim(c(-3.25,5))+
  geom_text_repel(aes(label=gene),size=2,fontface='italic',alpha=1,max.overlaps = 10,segment.size = 0.125)


g.cor.bab<-ggplot(imp.gene.lfc,aes(x=bat,y=baboon))+cust.theme()+xlab("bat LogFC")+ylab("baboon LogFC")+
  geom_point(data=imp.gene.lfc[!imp.gene.lfc$bat.baboon.sig=="FALSE FALSE",],aes(x=bat,y=baboon,colour=bat.baboon.sig),alpha=0.85,size=1.5)+
  geom_hline(yintercept=0,linetype='dotted')+geom_vline(xintercept=0,linetype='dotted')+
  theme(legend.position=c(4,2.5))+scale_colour_manual(values=c("#0C7C59","#444D5D","#EF6F6C"))#+#scale_color_manual(values=c("#888","#555"))+

# summary(imp.gene.lfc[imp.gene.lfc$gene %in% c(bat.sig1,pig.sig1),]$bat*imp.gene.lfc[imp.gene.lfc$gene %in% c(bat.sig1,pig.sig1),]$pig>0)#105/(105+45)
# summary(imp.gene.lfc[imp.gene.lfc$gene %in% c(bat.sig1,mac.sig1),]$bat*imp.gene.lfc[imp.gene.lfc$gene %in% c(bat.sig1,mac.sig1),]$mac>0)#105/(105+45)
# summary(imp.gene.lfc[imp.gene.lfc$gene %in% c(bat.sig1,bab.sig1),]$bat*imp.gene.lfc[imp.gene.lfc$gene %in% c(bat.sig1,bab.sig1),]$bab>0)#151/(151+77)

# spp.lps.heatmap /
#  (g.cor.pig+g.cor.mac+g.cor.bab)+ plot_layout(heights=c(1,2))


## upset plot
bat.sig.2<-imp.gene.lfc[imp.gene.lfc$batgene %in% bat.sig1 & imp.gene.lfc$bat>0,]$gene
bab.sig.2<-imp.gene.lfc[imp.gene.lfc$babgene %in% bab.sig1 & imp.gene.lfc$baboon>0,]$gene
mac.sig.2<-imp.gene.lfc[imp.gene.lfc$macgene %in% mac.sig1 & imp.gene.lfc$macaque>0,]$gene
pig.sig.2<-imp.gene.lfc[imp.gene.lfc$piggene %in% pig.sig1 & imp.gene.lfc$pig>0,]$gene

## ggplot2 compatible
gg.up.df<-data.frame(bat=rownames(imp.gene.lfc) %in% bat.sig.2,
                     baboon=rownames(imp.gene.lfc) %in% bab.sig.2,
                     macaque=rownames(imp.gene.lfc) %in% mac.sig.2,
                     pig=rownames(imp.gene.lfc) %in% pig.sig.2,row.names = rownames(imp.gene.lfc))


tidy_up_spp <- gg.up.df %>%
  as_tibble(rownames = "gene") %>%
  gather(spp, Member, -gene) %>%
  filter(Member) %>%
  dplyr::select(- Member)

tidy_up_spp<-tidy_up_spp %>%
  group_by(gene) %>%
  summarize(spp = list(spp))

tidy_up_spp$unique=TRUE
tidy_up_spp[lengths(tidy_up_spp$spp)>1,]$unique<-FALSE

upset.plot<-ggplot(tidy_up_spp,aes(x = spp,fill=unique)) + cust.theme()+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-1,size=2.75) +
  geom_bar(colour="#555555") +scale_fill_manual(values=c("#ccc","#888"))+
  scale_x_upset()+theme(panel.grid=element_blank(),legend.position=c(0.8,0.8))+
  scale_y_continuous(breaks = NULL, lim = c(0, 350), name = "")


## plot fig 1
(p7.pca+labs(tag="A")+g.volcano+labs(tag="B")) /
  (go.all.th.up+labs(tag="C")+dotplot.kegg+labs(tag="D")) /
  (as.ggplot(spp.lps.heatmap)+labs(tag="E")) /
  (upset.plot+labs(tag="F")+g.cor.mac+labs(tag="G"))+ 
  plot_layout(heights=c(1.4,1,1,1.4))

ggsave("fig1_new.svg",height=12,width=12)


## test overrepresentation of shared/unique LPS-responsive genes
gg.up.df2<-data.frame(bat=rownames(imp.gene.lfc) %in% bat.sig.2,
                      nonbat=rownames(imp.gene.lfc) %in% c(bab.sig.2,mac.sig.2,pig.sig.2))
table(gg.up.df2$bat,gg.up.df2$nonbat)
fisher.test(table(gg.up.df2$bat,gg.up.df2$nonbat))

### baboon comparison
# gg.up.df2<-data.frame(baboon=rownames(imp.gene.lfc) %in% bab.sig.2,
#                       nonbab=rownames(imp.gene.lfc) %in% c(bat.sig.2,mac.sig.2,pig.sig.2))
# fisher.test(table(gg.up.df2$baboon,gg.up.df2$nonbab))

bats.other.lps<-imp.gene.lfc[imp.gene.lfc$gene %in% bat.sig.2 & 
                               (imp.gene.lfc$gene %in% c(pig.sig.2,bab.sig.2,mac.sig.2)),]$gene

bats.only.lps<-imp.gene.lfc[imp.gene.lfc$gene %in% bat.sig.2 & 
                              !(imp.gene.lfc$gene %in% c(pig.sig.2,bab.sig.2,mac.sig.2)),]$gene

go.OR.batothersig<-clusterProfiler::simplify(enrichGO(gene = bats.other.lps,
                                                      universe = imp.gene.lfc$gene,#list of all genes
                                                      keyType = "SYMBOL",
                                                      OrgDb = organism,
                                                      ont = "BP",
                                                      pAdjustMethod = "BH",
                                                      pvalueCutoff = 0.05,
                                                      readable = TRUE),
                                             cutoff=0.7,by = "p.adjust",select_fun = min,
                                             measure = "Wang",semData = NULL)@result

go.OR.batonlysig<-clusterProfiler::simplify(enrichGO(gene = bats.only.lps,
                                                     universe =imp.gene.lfc$gene,#list of all genes
                                                     keyType = "SYMBOL",
                                                     OrgDb = organism,
                                                     ont = "BP",
                                                     pAdjustMethod = "BH",
                                                     pvalueCutoff = 0.05,
                                                     readable = TRUE),
                                            cutoff=0.7,by = "p.adjust",select_fun = min,
                                            measure = "Wang",semData = NULL)@result

plot_GO_fun(go.OR.batothersig)+labs(title="Bat and >= one other mammal")
ggsave("bat_other_GOoverrep.png",dpi=600,height=4,width=7)


#####################
# 05. SEX-SPECIFIC LPS EFFECTS

dds.trt.m <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.p7[,samples.p7$Sex=="M"]),
  colData = samples.p7[samples.p7$Sex=="M",],
  design= ~ Trtmt+Band))
summary((results(dds.trt.m, name="Trtmt_TH_vs_C0",alpha=alpha)))
res.th.vs.c0.m<- add_gene_info(data.frame(results(dds.trt.m, name="Trtmt_TH_vs_C0",alpha=alpha)))

dds.trt.f <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.p7[,samples.p7$Sex=="F"]),
  colData = samples.p7[samples.p7$Sex=="F",],
  design= ~ Trtmt+Band))
summary((results(dds.trt.f, name="Trtmt_TH_vs_C0",alpha=alpha)))
res.th.vs.c0.f<- add_gene_info(data.frame(results(dds.trt.f, name="Trtmt_TH_vs_C0",alpha=alpha)))

write.csv(res.th.vs.c0.m,"./phase7_paired_th_vs_c0_male_only.csv",row.names = FALSE)
write.csv(res.th.vs.c0.m,"./phase7_paired_th_vs_c0_female_only.csv",row.names = FALSE)

table(res.th.vs.c0.m$sig,res.th.vs.c0.f$sig)

## compare logFC
fm_slopes<-merge(res.th.vs.c0.f,res.th.vs.c0.m,by='gene')
fm_slopes$Trt_female=fm_slopes$log2FoldChange.x
fm_slopes$Trt_male=fm_slopes$log2FoldChange.y
fm_slopes$P_female=fm_slopes$padj.x
fm_slopes$P_male=fm_slopes$padj.y
fm_slopes$geneset<-"All genes"

### test correlation across all genes
cor.test(fm_slopes$Trt_female,fm_slopes$Trt_male)
### test correlation across LPS-responsive genes
cor.test(fm_slopes[fm_slopes$gene %in% res.th.vs.c0[res.th.vs.c0$sig,]$gene,]$Trt_female,
         fm_slopes[fm_slopes$gene %in% res.th.vs.c0[res.th.vs.c0$sig,]$gene,]$Trt_male)

g.fm_slopes<-ggplot(fm_slopes[fm_slopes$gene %in% res.th.vs.c0[res.th.vs.c0$sig,]$gene,],aes(x=Trt_female,y=Trt_male))+
  geom_point(data=fm_slopes,colour='#aaaaaa',size=1,aes(x=Trt_female,y=Trt_male))+
  geom_point(size=1,colour="#ed7272")+
  #geom_smooth(method='lm')+
  geom_abline(intercept=0,slope=1,colour='black',linewidth=0.5,linetype='dashed')+
  annotate(geom = 'text',x=8,y=-3,label="r = 0.961",colour="#ed7272")+
  annotate(geom = 'text',x=8,y=-4,label="r = 0.832",colour="#aaaaaa")+
  cust.theme()+xlab("Trt log2FC in females")+ylab("Trt log2FC in males")

### plot overrep
dotplot.sexdiff.fup<-plot_GO_overrep_upreg(res.th.vs.c0.f)
dotplot.sexdiff.mup<-plot_GO_overrep_upreg(res.th.vs.c0.m)
dotplot.sexdiff.fup+dotplot.sexdiff.mup

### now keep only concordant, sig genes
fm_slopes<-fm_slopes[fm_slopes$Trt_female*fm_slopes$Trt_male>0 & 
                       fm_slopes$gene %in% res.th.vs.c0[res.th.vs.c0$sig,]$gene,]

### test difference across LPS-responsive, concordant genes
t.test(fm_slopes[fm_slopes$gene %in% res.th.vs.c0[res.th.vs.c0$sig,]$gene,]$Trt_female,
       fm_slopes[fm_slopes$gene %in% res.th.vs.c0[res.th.vs.c0$sig,]$gene,]$Trt_male,paired=TRUE)

### plot differences in slope
fm_slopes$diff<-abs(fm_slopes$Trt_male)-abs(fm_slopes$Trt_female)
summary(fm_slopes$diff)

library(ggridges)
g.slope<-ggplot(fm_slopes,aes(x=diff,y=1))+
  geom_density_ridges2(scale = 0.85,alpha=1,colour='white',linewidth=0)+
  geom_boxplot(width=0.2,outlier.shape = NA,alpha=0.5,fill='white')+
  geom_vline(xintercept=0,linetype='dashed')+
  cust.theme()+theme(axis.text.y=element_blank(),legend.position='left',
                     axis.title.y=element_blank(),axis.ticks.y=element_blank(),legend.title=element_blank())+
  xlab('male abs. logFC - female abs. logFC')+
  coord_cartesian(xlim=c(-0.75,0.75))

ggsave('sex-specific_Trt_effect.png',
       plot=(g.fm_slopes+labs(tag="A",title='Sex correlation')+g.slope+labs(tag="B",title='Sex differences in LPS response'))+
         plot_layout(widths=c(1,1.2)),
       dpi=600,height=3.5,width=7)


#####################
# 06. VARIATION IN UNTREATED SAMPLES - BOTH SEXES

#test sex and sex effects across all C0 samples
dds.all.c0 <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.c0),
  colData = samples.all.c0,
  design= ~ Sex+Age+#lnNLR+
    Phase))
resultsNames(dds.all.c0)
summary((results(dds.all.c0, name="Age",alpha=alpha)))
summary((results(dds.all.c0, name="Sex_M_vs_F",alpha=alpha)))
# summary((results(dds.all.c0, name="lnNLR",alpha=alpha)))
# summary((results(dds.all.c0, name="SexM.Age",alpha=alpha)))
dds.trt.c0.age<-add_gene_info(data.frame(results(dds.all.c0, name="Age",alpha=alpha)))
dds.trt.c0.sex<-add_gene_info(data.frame(results(dds.all.c0, name="Sex_M_vs_F",alpha=alpha)))
write.csv(dds.trt.c0.age,"./results_files/phases3-7_C0_age_effect.csv",row.names = FALSE)
write.csv(dds.trt.c0.sex,"./results_files/phases3-7_C0_sex_effect.csv",row.names = FALSE)


#####################
# 07. VARIATION IN LPS-TREATED SAMPLES - BOTH SEXES

#test sex and sex effects across all TH samples
dds.all.th <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.th),
  colData = samples.all.th,
  design= ~ Sex+Age+#lnNLR+#including lnNLR or not has massive effect on sex-biased genes
    Phase))
resultsNames(dds.all.th)
summary((results(dds.all.th, name="Age",alpha=alpha)))
summary((results(dds.all.th, name="Sex_M_vs_F",alpha=alpha)))
#summary((results(dds.all.th, name="SexM.Age",alpha=alpha)))
# summary((results(dds.all.th, name="lnNLR",alpha=alpha)))

dds.trt.th.age<-add_gene_info(data.frame(results(dds.all.th, name="Age",alpha=alpha)))
dds.trt.th.sex<-add_gene_info(data.frame(results(dds.all.th, name="Sex_M_vs_F",alpha=alpha)))
write.csv(dds.trt.th.age,"./phases3-7_TH_age_effect.csv",row.names = FALSE)
write.csv(dds.trt.th.sex,"./phases3-7_TH_sex_effect.csv",row.names = FALSE)

## go plots for sex
go.th.sex.m<-plot_GO_overrep(dds.trt.th.sex,'up')+ labs(title="Sex - M up (logFC > 0)")
go.th.sex.f<-plot_GO_overrep(dds.trt.th.sex,'dn')+ labs(title="Sex - F up (logFC < 0)")

## go plots for age
go.th.age.up<-plot_GO_overrep(dds.trt.th.age,'up')+ labs(title="Age - up (logFC > 0)")
go.th.age.dn<-plot_GO_overrep(dds.trt.th.age,'dn')+ labs(title="Age - down (logFC < 0)")

## volcano pots
sex.vol<-plot_volcano(dds.trt.th.sex[-log10(dds.trt.th.sex$padj)<30,])+ labs(title="Sex-biased",tag="A")
age.vol<-plot_volcano(dds.trt.th.age) + labs(title="Age-associated",tag="B")

sex.vol + go.th.sex.m + go.th.sex.f + age.vol + go.th.age.up + go.th.age.dn+
  plot_layout(widths = c(1.5, 1),
              guides = "collect",
              design = "
              12
              13
              45
              46
              ")

ggsave('~/Documents/UMD_new/nih_rna/writeup/scripts_results_files/interindividual_var_sex_age.svg',dpi=600,height=8.5,width=9.5)

#are sex-biased genes overrepresented for LPS-affected genes
dds.trt.th.sex<-dds.trt.th.sex[!is.na(dds.trt.th.sex$padj),]
dds.trt.th.sex$LPS_affected<-"Not LPS-influenced"
dds.trt.th.sex[dds.trt.th.sex$gene %in% lps.genes,]$LPS_affected<-"LPS-influenced"
anno.set<-dds.trt.th.sex[abs(dds.trt.th.sex$log2FoldChange)>quantile(abs(dds.trt.th.sex$log2FoldChange),0.95) | 
                           -log10(dds.trt.th.sex$pvalue)>quantile(-log10(dds.trt.th.sex$pvalue),0.95),]
anno.set<-anno.set[-grep("^LOC",anno.set$gene),]
table(dds.trt.th.sex$sig,dds.trt.th.sex$LPS_affected)
fisher.test(table(dds.trt.th.sex$sig,dds.trt.th.sex$LPS_affected))



#####################
# 08. SEX-SPECIFIC AGE EFFECTS

## Perform sex-specific z-scaling
samples.all.th$ss.age<-samples.all.th$Est.Age
samples.all.th[samples.all.th$Sex=="M",]$ss.age<-scale(samples.all.th[samples.all.th$Sex=="M",]$Est.Age)
samples.all.th[samples.all.th$Sex=="F",]$ss.age<-scale(samples.all.th[samples.all.th$Sex=="F",]$Est.Age)

## run sex-specific models, WITH sex-specific z-scaling of ages
dds.all.th.f <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.th[,samples.all.th$Sex=="F"]),
  colData = samples.all.th[samples.all.th$Sex=="F",],
  design= ~ ss.age+Phase))
dds.all.th.m <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.th[,samples.all.th$Sex=="M",]),
  colData = samples.all.th[samples.all.th$Sex=="M",],
  design= ~ ss.age+Phase))
resultsNames(dds.all.th.m)
summary(results(dds.all.th.f, name="ss.age",alpha=alpha))
summary(results(dds.all.th.m, name="ss.age",alpha=alpha))
#summary(results(dds.all.th.m, name="STATUS_H_vs_B",alpha=alpha))

dds.trt.th.age.f<-add_gene_info(data.frame(results(dds.all.th.f, name="ss.age",alpha=alpha)))
dds.trt.th.age.m<-add_gene_info(data.frame(results(dds.all.th.m, name="ss.age",alpha=alpha)))

#plot_GO_overrep_upreg(dds.trt.th.age.f)+plot_GO_overrep_upreg(dds.trt.th.age.m)
#plot_GO_overrep_downreg(dds.trt.th.age.f)+plot_GO_overrep_downreg(dds.trt.th.age.m)
#write.csv(dds.trt.th.sex.f,"./results_files/phases3-7_TH_female_age_effect.csv",row.names = FALSE)
#write.csv(dds.trt.th.sex.m,"./results_files/phases3-7_TH_male_age_effect.csv",row.names = FALSE)

dds.th.age.ss<-merge(dds.trt.th.age.f,dds.trt.th.age.m,by='gene')
dds.th.age.ss.concord<-dds.th.age.ss[dds.th.age.ss$log2FoldChange.x*dds.th.age.ss$log2FoldChange.y>0 & 
                                       dds.th.age.ss$sig.x & dds.th.age.ss$sig.y,]
t.test(abs(dds.th.age.ss.concord$log2FoldChange.x),abs(dds.th.age.ss.concord$log2FoldChange.y),paired=TRUE)
(mean(2^abs(dds.th.age.ss.concord$log2FoldChange.y)-2^abs(dds.th.age.ss.concord$log2FoldChange.x)))
#sd(2^abs(dds.th.age.ss.concord$log2FoldChange.y)-2^abs(dds.th.age.ss.concord$log2FoldChange.x))/sqrt(120)


## run sex-specific models, WITHOUT sex-specific z-scaling of ages
dds.all.th.f <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.th[,samples.all.th$Sex=="F"]),
  colData = samples.all.th[samples.all.th$Sex=="F",],
  design= ~ Age+Phase))
dds.all.th.m <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.th[,samples.all.th$Sex=="M",]),
  colData = samples.all.th[samples.all.th$Sex=="M",],
  design= ~ Age+Phase))
resultsNames(dds.all.th.m)
summary(results(dds.all.th.f, name="Age",alpha=alpha))
summary(results(dds.all.th.m, name="Age",alpha=alpha))

dds.trt.th.age.f<-add_gene_info(data.frame(results(dds.all.th.f, name="Age",alpha=alpha)))
dds.trt.th.age.m<-add_gene_info(data.frame(results(dds.all.th.m, name="Age",alpha=alpha)))

dds.th.age<-merge(dds.trt.th.age.f,dds.trt.th.age.m,by='gene')
dds.th.age.concord<-dds.th.age[dds.th.age$log2FoldChange.x*dds.th.age$log2FoldChange.y>0 & 
                                 dds.th.age$sig.x & dds.th.age$sig.y,]
t.test(abs(dds.th.age.concord$log2FoldChange.x),abs(dds.th.age.concord$log2FoldChange.y),paired=TRUE)
(mean(2^abs(dds.th.age.concord$log2FoldChange.y)-2^abs(dds.th.age.concord$log2FoldChange.x)))
#sd(2^abs(dds.th.age.concord$log2FoldChange.y)-2^abs(dds.th.age.concord$log2FoldChange.x))/sqrt(120)

## run GO overrep. test of down-reg'd genes
go.concord.dn<-simplify(enrichGO(gene = dds.th.age.concord[dds.th.age.concord$log2FoldChange.x<0,]$gene,
                                 universe = dds.th.age$gene,#list of all genes
                                 keyType = "SYMBOL",
                                 OrgDb = organism,
                                 ont = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff = 1,
                                 readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                        measure = "Wang",semData = NULL)@result
plot_GO_fun(go.concord.dn)
ggsave("conc_downreg_MF_GOs.png",dpi=600,height=4,width=6)


## plot sex-specific age-related patterns, pre and post- z-scaling
dds.th.age.ss$sig<-factor(paste(dds.th.age.ss$sig.y,dds.th.age.ss$sig.x))
levels(dds.th.age.ss$sig)<-c(NA,"F","M","Both")

dds.th.age$sig<-factor(paste(dds.th.age$sig.y,dds.th.age$sig.x))
levels(dds.th.age$sig)<-c(NA,"F","M","Both")

g.age<-ggplot(samples.all.th,aes(x=Age,fill=Sex))+theme_minimal()+geom_density(alpha=0.5)+
  scale_fill_manual(values=SexPalette)+xlab("Age")+
  theme(axis.title.y=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),legend.position='none')
g.ss.age<-ggplot(samples.all.th,aes(x=ss.age,fill=Sex))+theme_minimal()+geom_density(alpha=0.5)+
  scale_fill_manual(values=SexPalette)+xlab("z-scaled age")+
  theme(axis.title.y=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),
        panel.grid=element_blank(),legend.position='none')

g.slope.th<-ggplot(dds.th.age,aes(x=log2FoldChange.x,y=log2FoldChange.y,colour=sig))+
  cust.theme()+geom_point(colour='grey',size=0.7)+
  geom_abline(intercept=0,slope=1,linewidth=0.25)+
  geom_point(data=dds.th.age[dds.th.age$sig.x | dds.th.age$sig.y,],
             size=1,alpha=1)+scale_colour_manual(values=c(SexPalette,"#b459ff"))+
  geom_point(data=dds.th.age[dds.th.age$sig.x & dds.th.age$sig.y,],#overlay genes sig in both sexes, so they're on top layer
             aes(x=log2FoldChange.x,y=log2FoldChange.y),colour='#b459ff',size=1.5)+
  theme(legend.position = 'right')+
  geom_hline(yintercept=0,linetype='dashed',linewidth=0.25)+
  geom_vline(xintercept=0,linetype='dashed',linewidth=0.25)+
  geom_smooth(data=dds.th.age[dds.th.age$sig.x & dds.th.age$sig.y,],
              aes(x=log2FoldChange.x,y=log2FoldChange.y),colour='black',method='lm',linewidth=1.5,se=FALSE)+#this just serves as a black outline
  geom_smooth(data=dds.th.age[dds.th.age$sig.x & dds.th.age$sig.y,],
              aes(x=log2FoldChange.x,y=log2FoldChange.y),colour='#b459ff',method='lm',linewidth=0.5,se=FALSE)+
  xlab("Age log2FC in LPS-treated females")+ylab("Age log2FC in LPS-treated males")+
  xlim(c(-4,2.5))+ylim(c(-4,2.5))

g.ss.slope.th<-ggplot(dds.th.age.ss,aes(x=log2FoldChange.x,y=log2FoldChange.y,colour=sig))+
  cust.theme()+geom_point(colour='grey',size=0.7)+
  geom_abline(intercept=0,slope=1,linewidth=0.25)+
  geom_point(data=dds.th.age.ss[dds.th.age$sig.x | dds.th.age.ss$sig.y,],
             size=1,alpha=1)+scale_colour_manual(values=c(SexPalette,"#b459ff"))+
  geom_point(data=dds.th.age.ss[dds.th.age.ss$sig.x & dds.th.age.ss$sig.y,],#overlay genes sig in both sexes, so they're on top layer
             aes(x=log2FoldChange.x,y=log2FoldChange.y),colour='#b459ff',size=1.5)+
  theme(legend.position = 'right')+
  geom_hline(yintercept=0,linetype='dashed',linewidth=0.25)+
  geom_vline(xintercept=0,linetype='dashed',linewidth=0.25)+
  geom_smooth(data=dds.th.age.ss[dds.th.age.ss$sig.x & dds.th.age.ss$sig.y,],
              aes(x=log2FoldChange.x,y=log2FoldChange.y),colour='black',method='lm',linewidth=1.5,se=FALSE)+#this just serves as a black outline
  geom_smooth(data=dds.th.age.ss[dds.th.age.ss$sig.x & dds.th.age.ss$sig.y,],
              aes(x=log2FoldChange.x,y=log2FoldChange.y),colour='#b459ff',method='lm',linewidth=0.5,se=FALSE)+
  xlab("Z-age log2FC in LPS-treated females")+ylab("Z-age log2FC in LPS-treated males")+
  xlim(c(-4,2.5))+ylim(c(-4,2.5))

g.age+labs(tag="A")+g.ss.age+labs(tag="B")+
  g.slope.th+theme(legend.position='none')+
  coord_cartesian(xlim=c(-3.5,2.5),ylim=c(-3.5,2.5))+
  g.ss.slope.th+coord_cartesian(xlim=c(-2.5,2),ylim=c(-2.5,2))+
  plot_layout(heights=c(0.35,1))

ggsave('sex_age.png',width=8,height=5)

