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
# 09. WGCNA module analyses of LPS-treated samples

library(sva)
library(WGCNA)
library(DESeq2)
library(flashClust)

# samples.all.th<-read.csv("~/Documents/UMD_new/nih_rna/writeup/samples_all_th.csv",row.names = 1)
# cts.all.th<-read.csv("~/Documents/UMD_new/nih_rna/writeup/cts_all_th.csv",row.names = 1)
# summary(colnames(cts.all.th)==rownames(samples.all.th))
# 
# dds.all.th <- DESeq(DESeqDataSetFromMatrix(
#   countData = round(cts.all.th),
#   colData = samples.all.th,
#   design= ~ Sex+Age+Phase))
# 
# vst<-data.frame(assay(vst(dds.all.th,blind=FALSE)))
# vst = ComBat(dat=vst, batch=samples.all.th$Phase, mod=NULL, par.prior=TRUE, prior.plots=FALSE)
# 
# write.csv(vst,"vst_phases3_7_TH_combat.csv",quote=FALSE)

samples.all.th<-read.csv("samples_all_th.csv",row.names = 1)
datExpr<-read.csv("vst_phases3_7_TH_combat.csv",row.names = 1)

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

#Create an object called "datTraits" that contains your trait data
datTraits = samples.all.th
head(datTraits)
table(rownames(datTraits)==rownames(datExpr)) #should return TRUE if datasets align correctly, otherwise your names are out of order

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
library(dplyr)
MEs0 <- orderMEs(MEs0)
module_order = names(MEs0) %>% gsub("ME","", .)

# Add treatment names
#MEs0$treatment = row.names(MEs0)
#MEs0$treatment = samples.all.th$Sex

datTraits[datTraits$Sex=="F",]$Sex<-2
datTraits[datTraits$Sex=="M",]$Sex<-1
datTraits<-datTraits[,c("Sex","Age")]
moduleTraitCor = cor(MEs0, datTraits, use= "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(datExpr))
moduleTraitPvalue[,1]=p.adjust(moduleTraitPvalue[,1])
moduleTraitPvalue[,2]=p.adjust(moduleTraitPvalue[,2])
moduleTraitPvalue<0.005
moduleTraitPvalue<0.05

textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitPvalue, 1), ")", sep= "")
dim(textMatrix)= dim(moduleTraitCor)
#par(mar= c(6, 8.5, 3, 3))

#display the corelation values with a heatmap plot
#INCLUE THE NEXT LINE TO SAVE TO FILE
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
#MEs0$Est.Age<-datTraits$ss.age
MEs0$Phase<-datTraits$Phase
MEs0$Band<-datTraits$Band

MEs0[MEs0$Sex==1,]$Sex="M"
MEs0[MEs0$Sex==2,]$Sex="F"

library(patchwork)
nrow(module_df[module_df$colors=="green",])
nrow(module_df[module_df$colors=="pink",])
nrow(module_df[module_df$colors=="yellow",])
nrow(module_df[module_df$colors=="blue",])

#run main scrpt first for pca plots

#MEs0<-MEs0[!duplicated(MEs0$Band),]
all.pca+stat_ellipse(level = 0.99)+theme(legend.position='left')+labs(tag="A")+g.pcdiff.th+theme(legend.position='none')+labs(tag="B")+
  ggplot(MEs0,aes(x=Est.Age,y=MEpink,colour=Sex))+cust.theme()+
  geom_point(alpha=0.75)+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+
  ylab("pink module (N = 165)")+labs(tag="C")+
  ggplot(MEs0,aes(x=Est.Age,y=MEgreen,colour=Sex))+cust.theme()+labs(tag="D")+
  geom_point(alpha=0.75)+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+
  ylab("green module (N = 424)")+theme(legend.position='none')+
  ggplot(MEs0,aes(x=Est.Age,y=MEyellow,colour=Sex))+cust.theme()+labs(tag="E")+
  geom_point(alpha=0.75)+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+
  ylab("yellow module (N = 732)")+theme(legend.position='none')+
  ggplot(MEs0,aes(x=Est.Age,y=MEblue,colour=Sex))+cust.theme()+labs(tag="F")+
  geom_point(alpha=0.75)+geom_smooth(method='lm')+scale_colour_manual(values=SexPalette)+
  ylab("blue module (N = 1653)")+theme(legend.position='none')+labs(tag="F")

ggsave("/Users/jackrayner/Documents/UMD_new/nih_rna/writeup/v2/PCA_wgcna_sexage.svg",dpi=600,height=5,width=9)
#ggsave("/Users/jackrayner/Documents/UMD_new/nih_rna/writeup/nodupes/PCA_wgcna_sexage.png",dpi=600,height=5,width=9)
#ggsave("/Users/jackrayner/Documents/UMD_new/nih_rna/writeup/v2/wgcna_sexage_sexspecific_zscale.svg",dpi=600,height=4.5,width=6)

g.go.green+g.go.yellow+g.go.blue


go.OR<-clusterProfiler::simplify(enrichGO(gene = module_df[module_df$colors=="green",]$gene_id,
                                          universe = module_df$gene_id,#list of all genes
                                          keyType = "SYMBOL",
                                          OrgDb = organism,
                                          ont = "BP",
                                          pAdjustMethod = "BH",
                                          pvalueCutoff = 0.05,
                                          readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                                 measure = "Wang",semData = NULL)@result
g.go.green<-plot_GO_fun(go.OR)#B cell stuff

go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="pink",]$gene_id,
                         universe = module_df$gene_id,#list of all genes
                         keyType = "SYMBOL",
                         OrgDb = organism,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                measure = "Wang",semData = NULL)@result
go.go.pink<-plot_GO_fun(go.OR)

go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="yellow",]$gene_id,
                         universe = module_df$gene_id,#list of all genes
                         keyType = "SYMBOL",
                         OrgDb = organism,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                measure = "Wang",semData = NULL)@result
g.go.yellow<-plot_GO_fun(go.OR)#T cell stuff

go.OR<-simplify(enrichGO(gene = module_df[module_df$colors=="blue",]$gene_id,
                         universe = module_df$gene_id,#list of all genes
                         keyType = "SYMBOL",
                         OrgDb = organism,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE),cutoff=0.7,by = "p.adjust",select_fun = min,
                measure = "Wang",semData = NULL)@result
g.go.blue<-plot_GO_fun(go.OR)#inflammation

go.go.pink+ggtitle("Pink module")+
  g.go.green+ggtitle("Green module")+
  g.go.blue+ggtitle("Blue module")+
  g.go.yellow+ggtitle("Yellow module")

ggsave("modules_GO.png",height=7,width=10)

library(lme4)
car::Anova(lmer(MEpink~Sex+Est.Age+(1|Phase),data=MEs0))
car::Anova(lmer(MEgreen~Sex*Est.Age+(1|Phase),data=MEs0),type="III")
car::Anova(lmer(MEyellow~Sex*Est.Age+(1|Phase),data=MEs0),type="III")
car::Anova(lmer(MEblue~Sex*Est.Age+(1|Phase),data=MEs0),type="III")

for (col in colnames(MEs0)){
  car::Anova(lm(col ~ Sex * Age,data=MEs0),type="III")
}

