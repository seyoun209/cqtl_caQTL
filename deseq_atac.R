library(data.table)
library(tidyverse)
library(DESeq2)
library(Rsamtools)
library(ggplot2)
library(magrittr)

atac_count <- fread("../peaks/CQTL_18_peakCounts.tsv",fill=TRUE, drop= c("V10"))
chr_nm <- sprintf("chr%s",seq(1:22))
atac_ct_subs <- atac_count[atac_count$chr %in% chr_nm,]

select(atac_ct_subs,contains("CTL"))

