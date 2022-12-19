library(ATACseqQC)
library(data.table)
library(dplyr)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(GenomicAlignments)


## prepare the example BAM files for importing
bamfiles_li <- list.files("../align",pattern = "nodups_filtered_sorted.bam$", full.names = TRUE)
bamFileLabels <- gsub("_18_nodups_filtered_sorted.bam", "", basename(bamfiles_li))

bamqc_fi <- NULL
for (i in 1:length(bamfiles_li)){
  print(i)
  bamFile <- bamfiles_li[i]
  bamqc_fi[[i]] <- bamQC(bamFile,outPath = NULL)
}
save(bamqc_fi,file="bamqc_fi")
