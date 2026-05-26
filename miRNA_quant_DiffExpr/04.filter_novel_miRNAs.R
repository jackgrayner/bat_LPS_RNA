nov<-read.csv("result_02_10_2025_t_20_06_39.tsv",h=T,sep="\t")
head(nov)
nrow(nov)

#get rid of low score, rfam alerts, lowly expressed miRNAs, and miRNAs with same seed as others
nov<-nov[nov$miRDeep2.score>10 & nov$rfam.alert=="-" & nov$total.read.count>10 & nov$example.miRBase.miRNA.with.the.same.seed=="-" & nov$significant.randfold.p.value=="yes",]
#assign random IDs
nov$Ph_ID<-paste0("pha-miR-",c(1000:(999+nrow(nov))))

nov$consensus.mature.sequence<-toupper(nov$consensus.mature.sequence)
nov$consensus.precursor.sequence<-toupper(nov$consensus.precursor.sequence)

nov.hairpin<-data.frame(nov[,c("Ph_ID","consensus.precursor.sequence")])
nov.hairpin$Ph_ID<-paste0(">",nov.hairpin$Ph_ID)
write.table(nov.hairpin,file="novel_pha_rnas_filtered_precursor.fa",sep="\n",row.names=FALSE,quote=FALSE,col.names=FALSE)

nov.mature<-data.frame(nov[,c("Ph_ID","consensus.mature.sequence")])
nov.mature$Ph_ID<-paste0(">",nov.mature$Ph_ID)
write.table(nov.mature,file="novel_pha_rnas_filtered_mature.fa",sep="\n",row.names=FALSE,quote=FALSE,col.names=FALSE)
