
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------------------------
# Read in eQTL results.
# ------------------------------------------------------------------------------

pvs0 = read.table("expression_vs_cn_with_pam50.txt", header=TRUE,
      sep="\t", as.is=TRUE)
dim(pvs0)
pvs0[1:2,]
summary(pvs0$p.value)

key0 = paste(pvs0$gene, pvs0$SNP, sep=":")

pvs = read.table("expression_vs_cn_with_pam50_qnorm.txt", header=TRUE,
          sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]
summary(pvs$p.value)

keys = paste(pvs$gene, pvs$SNP, sep=":")

length(intersect(key0, keys))
keyAll = union(key0, keys)
length(keyAll)

p0 = p1 = rep(5, length(keyAll))
p0[match(key0, keyAll)] = -log10(pvs0$p.value)
p1[match(keys, keyAll)] = -log10(pvs$p.value)

png("../figures2/expression_vs_cn_with_pam50_qnorm_vs_not.png",
      width=5, height=5, res=400, units="in")
par(bty="n", mar=c(5,4,1,1))
plot(p0, p1, xlab="-log10(p-value) Untransformed data",
             ylab="-log10(p-value) normalized data", pch=20,
             col=rgb(0.8,0.2,0.2,0.5))
abline(0, 1, col="grey", lwd=2)
dev.off()

table(pvs$p.value < 1e-20)
table(pvs$p.value < 1e-30)
table(pvs$p.value < 1e-40)
table(pvs$p.value < 1e-50)

# ------------------------------------------------------------------------------
# Read in location annotation.
# ------------------------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE, as.is=TRUE)
dim(infoE)
infoE[1:2,]

table(infoE$chr, useNA="ifany")

infoC = read.table("cn_info.txt", sep="\t", header=TRUE, as.is=TRUE)
dim(infoC)
infoC[1:2,]
table(infoC$chr, useNA="ifany")

# ------------------------------------------------------------------------------
# files
# ------------------------------------------------------------------------------

source("../../shared_code/eQTL_Plot.R")

inputs = "expression_vs_cn_qnorm.txt"
inputs = c(inputs, "expression_vs_cn_with_pam50_qnorm.txt")
inputs = c(inputs, "expression_vs_cn_with_ab_purity_pam50_qnorm.txt")

inputs

outputs = gsub(".txt", ".png", inputs, fixed=TRUE)

outputs

for(i in 1:3){
  
  inputi  = inputs[i]
  outputi = outputs[i]
  
  # ------------------------------------------------------------------------
  # read in eQTL results
  # ------------------------------------------------------------------------
  
  pvs = read.table(inputi, header=TRUE, sep="\t", as.is=TRUE)
  dim(pvs)
  pvs[1:2,]
  
  summary(pvs$p.value)
  
  table(pvs$p.value < 1e-30)
  table(pvs$p.value < 1e-40)
  table(pvs$p.value < 1e-50)
  table(pvs$p.value < 1e-60)
  table(pvs$p.value < 1e-80)
  
  # ------------------------------------------------------------------------
  # plot it
  # ------------------------------------------------------------------------
  
  table(pvs$gene %in% infoE$gene)
  table(pvs$SNP  %in% infoC$gene)
  
  ptplot   = pvs
  geneID   = match(ptplot$gene, infoE$gene)
  markerID = match(ptplot$SNP,  infoC$gene)
  scores   = ptplot$p.value
  scuts    = c(1e-20, 1e-30, 1e-40, 1e-60)
  cols     = c("green", "blue", "red", "black")
  eChr     = gsub("chr", "", infoE$chr)
  ePos     = 0.5*(infoE$start + infoE$end)
  mChr     = gsub("chr", "", infoC$chr)
  mPos     = 0.5*(infoC$start + infoC$end)
  chroms   = 1:22
  
  png(sprintf("../figures2/%s", outputi), width=7.5, height=9,
  res=200, units="in")
  
  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab="CN Location", ylab="Transcript Location",
  plot.hotspots=TRUE, hotspots.cut=10, score.type="p-value")
  dev.off()
  
}

q(save="no")

