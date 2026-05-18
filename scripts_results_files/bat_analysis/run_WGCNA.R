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

# nb. script assumes load_data_functions_summarise.R already run.

#####################
# 09. WGCNA module analyses of LPS-treated samples

library(sva)
library(WGCNA)
library(DESeq2)
library(flashClust)

summary(colnames(cts.all.th)==rownames(samples.all.th))

dds.all.th <- DESeq(DESeqDataSetFromMatrix(
  countData = round(cts.all.th),
  colData = samples.all.th,
  design= ~ Sex+Age+Phase))

#take variance stabilised counts
vst<-data.frame(assay(vst(dds.all.th,blind=FALSE)))
#renmove phase effects
vst = ComBat(dat=vst, batch=samples.all.th$Phase, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
write.csv(vst,"vst_phases3_7_TH_combat.csv",quote=FALSE)

datExpr<-read.csv("vst_phases3_7_TH_combat.csv",row.names = 1)
datTraits = samples.all.th
table(rownames(datTraits)==rownames(datExpr)) 

rowvars.dat<-(rowVars(as.matrix(datExpr)))
summary(rowvars.dat>quantile(rowvars.dat,0.5))
datExpr<-datExpr[rowvars.dat>quantile(rowvars.dat,0.5),]

datExpr = as.data.frame(t(datExpr)) 
dim(datExpr)
gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

#If the last statement returns TRUE, all genes have passed the cuts. If not, we remove the offending genes and samples from the data with the following:
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr)[!gsg$goodGenes], collapse= ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse=", ")))
  datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
}

gsg = goodSamplesGenes(datExpr, verbose = 3)
gsg$allOK

powers = c(c(1:10), seq(from =10, to=30, by=1)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr, powerVector=powers, verbose =5, networkType="signed") #call network topology analysis function
sizeGrWindow(9,5)
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")

enableWGCNAThreads()
picked_power = 6
temp_cor <- cor       
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(datExpr,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "signed",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 30,
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "ER",
                          
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)

mergedColors = labels2colors(netwk$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(
  netwk$dendrograms[[1]],
  mergedColors[netwk$blockGenes[[1]]],
  "Module colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05 )

module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

write_delim(module_df,
            file = "gene_modules.txt",
            delim = "\t")

# Get Module Eigengenes per cluster
MEs0 <- moduleEigengenes(datExpr, mergedColors)$eigengenes

# Reorder modules so similar modules are next to each other
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

datTraits[datTraits$Sex=="F",]$Sex<-2
datTraits[datTraits$Sex=="M",]$Sex<-1
datTraits<-datTraits[,c("Sex","Age")]
moduleTraitCor = cor(MEs0, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))
moduleTraitPvalue[,1]=p.adjust(moduleTraitPvalue[,1])
moduleTraitPvalue[,2]=p.adjust(moduleTraitPvalue[,2])
moduleTraitPvalue<0.05
moduleTraitPvalue.df<-data.frame(moduleTraitPvalue)
sig.mods<-rownames(
  moduleTraitPvalue.df[(moduleTraitPvalue.df$Sex<0.05 | moduleTraitPvalue.df$Age<0.05),]
)

textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)

#display the corelation values with a heatmap plot
pdf('heatmap.pdf')
labeledHeatmap(Matrix= moduleTraitCor,
               xLabels= names(datTraits),
               yLabels= names(MEs0),
               ySymbols= names(MEs0),
               colorLabels= FALSE,
               colors= blueWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim= c(-1,1),
               main= paste("Module-trait relationships"))
moduleTraitPvalue<-data.frame(moduleTraitPvalue)
dev.off()

MEs0$Sex<-datTraits$Sex
MEs0$Est.Age<-datTraits$Age

MEs0[MEs0$Sex==1,]$Sex="M"
MEs0[MEs0$Sex==2,]$Sex="F"

#run main scrpt first for pca plots
MEs0$Est.Age<-samples.all.th$Est.Age

sig.mods
nrow(module_df[module_df$colors==gsub("ME","",sig.mods[1]),])
nrow(module_df[module_df$colors==gsub("ME","",sig.mods[2]),])
nrow(module_df[module_df$colors==gsub("ME","",sig.mods[3]),])
nrow(module_df[module_df$colors==gsub("ME","",sig.mods[4]),])

all.pca+stat_ellipse(level = 0.99,aes(colour=trt),show.legend = FALSE,linewidth=0.4)+scale_colour_manual(values=TrtPalette)+
  theme(legend.position='left')+labs(tag="A",fill='Treatment',subtitle="PCA - all")+
  g.pcdiff.th+theme(legend.position='none')+labs(tag="B",x="Est.Age",y="PC1",subtitle = "PC1 - all")+
  p7.pcdiff+labs(tag="D",subtitle=expression(PC1[LPS] - PC1[untreated]))+
  ggplot(MEs0,aes(x=Est.Age,y=MEblue,fill=Sex,colour=Sex))+cust.theme()+theme(legend.position='none')+
  geom_point(alpha=0.75,shape=21,stroke=0.25,aes(fill=Sex),size=2,colour='black')+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+scale_fill_manual(values=SexPalette)+
  ylab("Blue module (N = 1648)")+labs(tag="C",subtitle="\"Inflammation\" module")+
  #ggplot(MEs0,aes(x=Est.Age,y=MEmagenta,fill=Sex,colour=Sex))+cust.theme()+labs(tag="D")+
  #geom_point(alpha=0.75,shape=21,stroke=0.25,aes(fill=Sex),size=2,colour='black')+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+scale_fill_manual(values=SexPalette)+
  #ylab("Magenta module (N = 97)")+theme(legend.position='none')+
  ggplot(MEs0,aes(x=Est.Age,y=MEblack,fill=Sex,colour=Sex))+cust.theme()+labs(tag="E",subtitle="\"B cell activation\" module")+
  geom_point(alpha=0.75,shape=21,stroke=0.25,aes(fill=Sex),size=2,colour='black')+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+scale_fill_manual(values=SexPalette)+
  ylab("Black module (N = 232)")+theme(legend.position='none')+
  ggplot(MEs0,aes(x=Est.Age,y=MEred,fill=Sex,colour=Sex))+cust.theme()+labs(tag="F",subtitle="\"T cell activation\" module")+
  geom_point(alpha=0.75,shape=21,stroke=0.25,aes(fill=Sex),size=2,colour='black')+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+scale_fill_manual(values=SexPalette)+
  ylab("Red module (N = 388)")+theme(legend.position='none')+labs(tag="F")

ggsave("PCA_wgcna_sexage.svg",dpi=600,height=5.5,width=9)
#ggsave("PCA_wgcna_sexage_nodupes.svg",dpi=600,height=5,width=9)


go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="blue",]$gene_id,
                                          universe = module_df$gene_id,#list of all genes
                                          keyType = "SYMBOL",
                                          OrgDb = organism,
                                          ont = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                                 measure = "Wang",semData = NULL)@result
g.go.blue<-plot_GO_fun(go.OR)#B cell stuff

go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="magenta",]$gene_id,
                         universe = module_df$gene_id,#list of all genes
                         keyType = "SYMBOL",
                         OrgDb = organism,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 1,
                         readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                measure = "Wang",semData = NULL)@result
g.go.magenta<-plot_GO_fun(go.OR)

go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="black",]$gene_id,
                         universe = module_df$gene_id,#list of all genes
                         keyType = "SYMBOL",
                         OrgDb = organism,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                measure = "Wang",semData = NULL)@result
g.go.black<-plot_GO_fun(go.OR)#T cell stuff

go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="red",]$gene_id,
                         universe = module_df$gene_id,#list of all genes
                         keyType = "SYMBOL",
                         OrgDb = organism,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                measure = "Wang",semData = NULL)@result
g.go.red<-plot_GO_fun(go.OR)#inflammation

g.go.blue+ggtitle("Blue module")+
  #g.go.magenta+ggtitle("Magenta module")+
  g.go.black+ggtitle("Black module")+
  g.go.red+ggtitle("Red module")+plot_layout(nrow=2)

ggsave("modules_GO.png",height=6,width=10)
ggsave("modules_GO.pdf",height=7,width=10)

library(lme4)
car::Anova(lm(MEblue~Sex*Est.Age,data=MEs0),type="II")
car::Anova(lm(MEmagenta~Sex*Est.Age,data=MEs0),type="II")
car::Anova(lm(MEblack~Sex*Est.Age,data=MEs0),type="II")
car::Anova(lm(MEred~Sex*Est.Age,data=MEs0),type="III")

car::Anova(lm(MEblue~Est.Age,data=MEs0[MEs0$Sex=="M",]),type="II")
car::Anova(lm(MEblue~Est.Age,data=MEs0[MEs0$Sex=="F",]),type="II")
car::Anova(lm(MEmagenta~Est.Age,data=MEs0[MEs0$Sex=="M",]),type="II")
car::Anova(lm(MEblack~Est.Age,data=MEs0[MEs0$Sex=="M",]),type="II")
car::Anova(lm(MEred~Est.Age,data=MEs0[MEs0$Sex=="M",]),type="II")



sig.mods1<-cbind(MEs0$MEblue,MEs0$MEmagenta,MEs0$MEblack,MEs0$MEred)
man1<-manova(sig.mods1~Sex*Est.Age,data=MEs0)
summary(man1)

for (col in colnames(MEs0)){
  car::Anova(lm(col ~ Sex * Age,data=MEs0),type="III")
}
