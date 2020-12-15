
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoE)
infoE[1:2,]

infoM = read.table("methylation_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoM)
infoM[1:2,]

eChr = gsub("chr", "", infoE$chr)
mChr = infoM$Chromosome

table(eChr, useNA="ifany")
table(mChr, useNA="ifany")

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

inputs = "expression_vs_methylation.txt"
inputs = c(inputs, "expression_vs_methylation_pam50.txt")
inputs = c(inputs, "expression_vs_methylation_absolute_purity.txt")
inputs = c(inputs, "expression_vs_methylation_ab_purity_pam50.txt")
inputs = c(inputs, "expression_vs_methylation_rm_cn.txt")
inputs = c(inputs, "expression_vs_methylation_rm_cn_pam50.txt")
inputs = c(inputs, "expression_vs_methylation_rm_cn_ab_purity.txt")
inputs = c(inputs, "expression_vs_methylation_rm_cn_ab_purity_pam50.txt")

labels = c("baseline", "pam50", "ab_purity", "ab_purity_pam50")
labels = c(labels, c("cn", "cn_pam50", "cn_ab_purity", "cn_ab_purity_pam50"))

pcuts = c(0, 10^(seq(-100, -40, by=10)))
np    = length(pcuts) - 1

tbl1 = tbl2 = matrix(nrow=length(inputs), ncol=np)

for(i in 1:length(inputs)){
  
  ff1 = inputs[i]
  
  cat(i, date(), ff1, "\n")
  
  pvs = read.table(ff1, header=TRUE, sep="\t", as.is=TRUE)
  dim(pvs)
  pvs[1:2,]

  summary(pvs$p.value)

  pvs = pvs[pvs$p.value < 1e-40,]
  dim(pvs)
  pvs[1:2,]

  peChr = eChr[match(pvs$gene, infoE$gene)]
  pmChr = mChr[match(pvs$SNP,  infoM$Composite.Element.REF)]

  for(j in 1:np){
    
    ww1 = which(pvs$p.value < pcuts[j+1] & pvs$p.value >= pcuts[j])
    tbl1[i,j] = length(ww1)
    tbl2[i,j] = length(which(peChr[ww1] == pmChr[ww1]))/length(ww1)
  }
}

pp1s = pcuts[1:(length(pcuts)-1)]
pp2s = pcuts[2:length(pcuts)]

rownames(tbl1) = rownames(tbl2) = labels
colnames(tbl1) = colnames(tbl2) = sprintf("[%.0e, %.0e)", pp1s, pp2s)

tbl1
tbl2

tbl1 = cbind(tbl1, rowSums(tbl1))
colnames(tbl1)[ncol(tbl1)] = "Total"

write.table(tbl1, sep=" & ", eol = "\\\\\n", quote=FALSE)

localPercent = rowSums(tbl1[,1:5]*tbl2[,1:5])/tbl1[,6]
tbl2 = cbind(tbl2, localPercent)
tbl2 = round(tbl2, 3)

write.table(tbl2, sep=" & ", eol = "\\\\\n", quote=FALSE)

sort(rowSums(tbl1))

write.table(tbl1, file = "eQTL_hits.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

write.table(tbl2, file = "eQTL_local_percentage.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)


q(save = "no")

