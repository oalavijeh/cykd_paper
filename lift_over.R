library(liftover)
library(gwascat)
library(data.table)
meta <- fread("Documents/Projects/pkd/meta_analysis/hetero_no_overlap_gel_finngenr8_ukbb_jbb_meta_analysis_for_manhattan_fuma.tbl")
meta$id <- paste0(meta$chr,"_",meta$pos,"_",meta$ref,"_",meta$alt)
chr <- meta$chr
start <- meta$pos
snp.name <- meta$id
chain.file.path <- "Documents/Reference/hg38ToHg19.over.chain"
example.38.gr<-GRanges(
  seqname=Rle(paste("chr",chr,sep="")),
  ranges=IRanges(start=start,end=start),
  snp.name=snp.name)

c<-import.chain(chain.file.path) ## e.g. hg38ToHg19.over.chain
example.37.gr<-unlist(liftOver(example.38.gr,c))  
example.37.gr <- as.data.frame(example.37.gr)
meta37 <- merge(x=meta[,c(1:14,20)], y=example.37.gr, by.x = "id", by.y = "snp.name", all.x=T, all.y=F) 
data.frame(colnames(meta37))
meta37 <- meta37[,c(1,2,17,4:15)]
write.table(meta37, "Documents/Projects/pkd/meta_analysis/b37_lifted_over_hetero_no_overlap_gel_finngenr8_ukbb_jbb_meta_analysis_cystic_for_manhattan_fuma.tbl", col.names = T, row.names = F, quote = F, sep = '\t')
