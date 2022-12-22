library(ATACseqQC)
library(data.table)
library(dplyr)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(GenomicAlignments)
library(ggpubr)


## prepare the example BAM files for importing
bamfiles_li <- list.files("../align",pattern = "nodups_filtered_sorted.bam$", full.names = TRUE)
bamFileLabels <- gsub("_18_nodups_filtered_sorted.bam", "", basename(bamfiles_li))
bamFileLabels <- gsub("CQTL_","",bamFileLabels)

bamqc_fi <- NULL
for (i in 1:length(bamfiles_li)){
  print(i)
  bamFile <- bamfiles_li[i]
  bamqc_fi[[i]] <- bamQC(bamFile,outPath = NULL)
}
#save(bamqc_fi,file="bamqc_fi")
load("bamqc_fi")

totalQNAMEs <- sapply(bamqc_fi,'[[',1)
duplicateRate <- sapply(bamqc_fi,'[[',2)
properPairRate <- sapply(bamqc_fi,'[[',4)
unmappedRate <- sapply(bamqc_fi,'[[',5)
nonRedundantFraction <- sapply(bamqc_fi,'[[',8)

qc_subset <- cbind(totalQNAMEs,duplicateRate,properPairRate,unmappedRate,nonRedundantFraction) %>% as.data.frame()
qc_subset$samplenm <- bamFileLabels

totalread <- ggplot(data=qc_subset,aes(x=samplenm,y=totalQNAMEs)) +
  geom_bar(stat="identity",fill="#E64B35B2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=totalQNAMEs), vjust=1.6, color="black", size=3.5)

dups <- ggplot(data=qc_subset,aes(x=samplenm,y=duplicateRate)) +
  geom_bar(stat="identity",fill="#4DBBD5B2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=round(duplicateRate, digits = 4)), vjust=1.6, color="black", size=3.5)

pairrate <-ggplot(data=qc_subset,aes(x=samplenm,y=properPairRate)) +
  geom_bar(stat="identity",fill="#8491B4B2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=round(properPairRate, digits = 4)), vjust=1.6, color="black", size=3.5)

unmappedrate <- ggplot(data=qc_subset,aes(x=samplenm,y=unmappedRate)) +
  geom_bar(stat="identity",fill="#00A087B2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=round(unmappedRate, digits = 4)), vjust=1.6, color="black", size=3.5)

nonreduc <- ggplot(data=qc_subset,aes(x=samplenm,y=nonRedundantFraction)) +
  geom_bar(stat="identity",fill="#7E6148B2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label=round(nonRedundantFraction, digits = 4)), vjust=1.6, color="black", size=3.5)

png(paste0("01.figure_qc/qc_subset",".png"))
ggarrange(totalread,dups,pairrate,unmappedrate,nonreduc,ncol = 2, nrow = 3)
dev.off()

# fragsize <- NULL
# for (i in 1:length(bamfiles_li)) {
#   print(i)
#   bamFile <- bamfiles_li[i]
#   fragsize[[i]] <- fragSizeDist(bamFile, bamFileLabels[i])
#   png(paste0("01.figure_qc/",bamFileLabels[i],".png"))
#   fragSizeDist(bamFile, bamFileLabels[i])
#   dev.off()
# }

lib_complexity <- NULL
for (i in 2:length(bamfiles_li)) {
  print(i)
  bamFile <- bamfiles_li[i]
  lib_complexity[[i]] <- estimateLibComplexity(readsDupFreq(bamFile))
  png(paste0("01.figure_qc/",bamFileLabels[i],"libcomplexity",".png"))
  estimateLibComplexity(readsDupFreq(bamFile))
  dev.off()
}
save(lib_complexity,file=lib_complexity)

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), fix = "start", 1)
TSSs
chr_nm <- sprintf("chr%s",c(seq(1:22),"X","Y"))
gals <- lapply(bamfiles_li,function(bamfile){
  readBamFile(bamFile=bamfile, tag=character(0),
              which=GRanges(chr_nm, IRanges(1, 1e6)),
              asMates=FALSE)
})
bamFileLabels <- gsub("CQTL_","",bamFileLabels)
names(gals) <- bamFileLabels

txs <- transcripts(TxDb.Hsapiens.UCSC.hg38.knownGene)
png(paste0("01.figure_qc/","similarity_n3",".png"))
plotCorrelation(GAlignmentsList(gals), txs,seqlev=chr_nm,type = c("heatmap"),cexRow=0.8,cexCol=0.8)
dev.off()

png(paste0("01.figure_qc/","similarity_n3_pca",".png"))
plotCorrelation(GAlignmentsList(gals), txs,seqlev=chr_nm,type = c("PCA"))
dev.off()


