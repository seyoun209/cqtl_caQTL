library(ChIPseeker)
library(ChIPQC)
library(rtracklayer)
library(DT)
library(dplyr)
library(tidyr)
library(soGGi)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)



blkList <- import.bed("ATAC_Data/ATAC_blacklists/ENCFF001TDO.bed.gz")
peaks <- dir("../peaks", pattern = "*.narrowPeak", 
             full.names = TRUE)
myPeaks <- lapply(peaks, ChIPQC:::GetGRanges, simple = TRUE)
do.call(rbind,strsplit(do.call(rbind,strsplit(peaks,"/"))[,3],"_18"))[,1]
names(myPeaks) <- c(do.call(rbind,strsplit(do.call(rbind,strsplit(peaks,"/"))[,3],"_18"))[,1])

filteredAnno <- NULL
for (i in 1:length(peaks)){
  filteredAnno[[i]] <- annotatePeak(peaks[i], TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene)
}

for (i in 1:length(peaks)){
  png(paste0("01.figure_qc/anno_pichart_",names(filteredAnno)[i],".png"))
  plotAnnoPie(filteredAnno[[i]])
  text(0.5,0.7,names(filteredAnno)[i],font=2)
  dev.off()
}

png(paste0("01.figure_qc/anno_bar_all",".png"))
plotAnnoBar(filteredAnno)
dev.off()


### Diff atac-seq
myGRangesList<-GRangesList(myPeaks)
Group <- factor(c(do.call(rbind,strsplit(do.call(rbind,strsplit(peaks,"/"))[,3],"_"))[,3]))
consensusToCount <- soGGi:::runConsensusRegions(GRangesList(myPeaks), "none")

