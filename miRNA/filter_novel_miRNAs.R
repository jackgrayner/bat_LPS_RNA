nov<-read.csv("~/Documents/UMD_new/nih_rna/miRNA/novel_pha_rnas.csv",h=T)
head(nov)
nrow(nov)

#get rid of low score, rfam alerts, lowly expressed miRNAs, and miRNAs with same seed as others
nov<-nov[nov$miRDeep2.score>10 & nov$rfam.alert=="-" & nov$total.read.count>10 & nov$example.miRBase.miRNA.with.the.same.seed=="-",]
#assign random IDs
nov$Ph_ID<-paste0("phha-miR-",c(1:nrow(nov)))

nov$consensus.mature.sequence<-toupper(nov$consensus.mature.sequence)
nov$consensus.precursor.sequence<-toupper(nov$consensus.precursor.sequence)
write.csv(nov[,c("Ph_ID","consensus.precursor.sequence")],file="~/Documents/UMD_new/nih_rna/miRNA/novel_pha_rnas_filtered_precursor.csv",row.names=FALSE,quote=FALSE)
write.csv(nov[,c("Ph_ID","consensus.mature.sequence")],file="~/Documents/UMD_new/nih_rna/miRNA/novel_pha_rnas_filtered_mature.csv",row.names=FALSE,quote=FALSE)
