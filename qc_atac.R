library(ATACseqQC)
library(data.table)
library(dplyr)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)


## prepare the example BAM files for importing
bamfiles_li <- list.files("../align",pattern = "nodups_filtered_sorted.bam$", full.names = TRUE)
bamFileLabels <- gsub("_18_nodups_filtered_sorted.bam", "", basename(bamfiles_li))

bamqc_fi <- NULL
for (i in 1:length(bamfiles_li)){
  print(i)
  bamFile <- bamfiles_li[i]
  bamqc_fi[[i]] <- bamQC(bamFile)
}
save(bamqc_fi,file="bamqc_result")

fragsize <- NULL
for (i in 1:length(bamfiles_li)) {
  print(i)
  bamFile <- bamfiles_li[i]
  fragsize[[i]] <- fragSizeDist(bamFile, bamFileLabels[i])
  png(paste0("01.figure_qc/",bamFileLabels[i],".png"))
  fragSizeDist(bamFile, bamFileLabels[i])
  dev.off()
}

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
TSSs <- resize(genes(TxDb.Hsapiens.UCSC.hg38.knownGene), fix = "start", 1)
TSSs

gals <- lapply(bamfiles_li, function(bamfile){
  readBamFile(bamFile=bamfile, tag=character(0),
              which=GRanges("chr1", IRanges(1, 1e6)),
              asMates=FALSE)
})
