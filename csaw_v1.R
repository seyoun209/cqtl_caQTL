library(csaw)
library(BSgenome.Hsapiens.UCSC.hg38)
bamFile <- list.files("../align",pattern = "nodups_filtered_sorted.bam$", full.names = TRUE)
bamFileLabels <- gsub("_18_nodups_filtered_sorted.bam", "", basename(bamFile))

frag.lens=c();
for (i in 1:length(bamFile)){
  print(i)
  out = getPESizes(bamFile[i]);
  frag.lens=c(frag.lens,mean(out$sizes));
  frag.sizes <- out$sizes[out$sizes<=2000];
  hist(frag.sizes, breaks=50, xlab="Fragment sizes (bp)",ylab="Frequency", main="", col="grey80");
  abline(v=1500, col="red");
}
dev.off();
#save(frag.lens, file="frag.lens.Rdata");
load("frag.lens.Rdata");
pe.param = readParam(max.frag=1500, pe="both",minq=20);
#framgment length
multi.frag.len=list(frag.lens,NA);
#binding site length
win.width=10;
#spacing should notbe larger than ext/2 for analyses with small windows. 
#If ext is also very small, spacing should be set to width to avoid loading too many small windows.
data = windowCounts(bamFile, ext=multi.frag.len, filter=700, spacing=win.width, param=pe.param, width=win.width);
merged <- mergeWindows(rowRanges(data), tol=100L);
my.regions = merged$region;
save(my.regions,file="regions.Rdata");

##change chromosome names
my.regions=GRanges(seqnames(my.regions),IRanges(start(my.regions),end(my.regions)));

##counts
reg.counts <- regionCounts(bamFile, my.regions, ext=multi.frag.len, param=pe.param);
counts=assay(reg.counts);
save(counts,file="CSAW_100L_counts.Rdata");
