# Permutation testing using case-control status

library(data.table)
library(coin)

sv="INV"
svlen=0.05
pheno ="puv"

# Read in SV summary 
both<- read.delim(paste0("Rare_StructuralVariants_puv_ancestry_matched_controls_",bed,"_24814x_participants.txt"), header=T, sep="\t")
both <- both[!(both$CHROM=="chrX" | both$CHROM=="chrY"),]
both <- both[,c(8,9,10,4,1,2,3,7,11,5)]
both$SVlengthkb <-an(both$LENGTH.Mb.)*1000
both <- both[,-c(9)]
both <- na.omit(both)
both <- both[!duplicated(both[,c(1,4,5)]),]

# Read in ped file
ped <- fread("~/re_gecip/renal/mchan/PUV/pca/puv_ancestry_matched_controls_unrelated_aggv2_final.ped", data.table=F)
ped <- na.omit(ped)

# Create data frame using presence of rare SV as outcome 
k <- ped[,c(1,6)]
l <-both[which(both$CONSEQUENCE==sv & both$SVlengthkb>svlen),]
ids <- unique(l$Part_ID)
k$outcome <- ifelse(k$FID %in% ids,1,0)

# Permutation testing
all <- independence_test(outcome~PHENO, data=k, alternative="two.sided")

# By cCRE
perm_test <- function(x, cre){
  m <- x[grep(paste0(cre),as.vector(x$GeneSymbol)),]
  cre_ids <- unique(m$Part_ID)
  k$outcome <- ifelse(k$FID %in% cre_ids,1,0)
  independence_test(outcome~PHENO, data=k, alternative="two.sided")
}

perm_test(l,"dELS")
perm_test(l,"pELS")
perm_test(l,"PLS")
perm_test(l,"CTCF-only")
perm_test(l,"DNase-H3K4me3")

