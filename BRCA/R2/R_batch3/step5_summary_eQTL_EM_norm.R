
# ------------------------------------------------------------------------------
# Read in eQTL results.
# ------------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

pvs0 = read.table("expression_vs_methylation_rm_cn_pam50.txt",
            header=TRUE, sep="\t", as.is=TRUE)
dim(pvs0)
pvs0[1:2,]
summary(pvs0$p.value)

key0 = paste(pvs0$gene, pvs0$SNP, sep=":")

pvs = read.table("expression_vs_methylation_rm_cn_pam50_qnorm.txt",
            header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

# show the relation between regression coefficient and p-value cutoff
summary(pvs$p.value)

png("../figures2/p-value_vs_beta.png", width=6, height=4,
    res=400, units="in")
par(bty="n", mar=c(5,4,1,1))
smoothScatter(pvs$beta, -log10(pvs$p.value), xlab="beta", ylab="-log10(p-value)")
dev.off()

summary(pvs$p.value)
table(pvs$p.value < 1e-60)
table(abs(pvs$beta) > 0.75)

# ------------------------------------------------------------
# read in gene annotation
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
                   as.is=TRUE)
dim(infoE)
infoE[1:2,]

# ------------------------------------------------------------------------
# select genes that appear at least 100 times in the beta > 0.75 list
# ------------------------------------------------------------------------

genes = pvs$gene[which(abs(pvs$beta) > 0.75)]

tg = table(genes)
summary(as.numeric(tg))

gene100 = names(tg)[which(as.numeric(tg) >= 100)]
length(gene100)

table(gene100 %in% infoE$gene)

gene100 = infoE[match(gene100, infoE$gene),]
dim(gene100)
gene100[1:5,]

table(gene100$chr)
table(is.na(gene100$ensembl))

cat(gene100$ensembl, sep="\n")

# ------------------------------------------------------------------------
# compare p-values with or without quantile normalization
# ------------------------------------------------------------------------

keys = paste(pvs$gene, pvs$SNP, sep=":")

length(intersect(key0, keys))
keyAll = union(key0, keys)
length(keyAll)

summary(pvs0$p.value)
summary(pvs$p.value)

p0 = p1 = rep(39, length(keyAll))
p0[match(key0, keyAll)] = -log10(pvs0$p.value)
p1[match(keys, keyAll)] = -log10(pvs$p.value)

png("../figures2/expression_vs_methylation_rm_cn_pam50_qnorm_vs_not.png",
    width=5, height=5,
res=400, units="in")
par(bty="n", mar=c(5,4,1,1))
plot(p0, p1, xlab="-log10(p-value) Untransformed data",
ylab="-log10(p-value) normalized data", pch=20,
col=rgb(0.8,0.2,0.2,0.5))
abline(0, 1, col="grey", lwd=2)
dev.off()

table(pvs$p.value < 1e-50)
table(pvs$p.value < 1e-60)

# --------------------------------------------------------------
# files
# --------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

inputs = "expression_vs_methylation_qnorm.txt"
inputs = c(inputs, "expression_vs_methylation_rm_cn_qnorm.txt")
inputs = c(inputs, "expression_vs_methylation_rm_cn_pam50_qnorm.txt")
inputs = c(inputs, "expression_vs_methylation_rm_cn_ab_purity_pam50_qnorm.txt")

inputs

outputs = gsub(".txt", ".png", inputs, fixed=TRUE)

outputs

# --------------------------------------------------------------
# read in location annotation
# --------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoE)
infoE[1:2,]

table(infoE$chr, useNA="ifany")

infoM = read.table("methylation_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoM)
infoM[1:2,]
table(infoM$Chromosome, useNA="ifany")

# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# --------------------------------------------------------------
# iteratively summarizing results
# --------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

source("~/research/TCGA/_Sun_MethyE/shared_code/eQTL_Plot.R")

for(i in 1:4){
  
  inputi  = inputs[i]
  outputi = outputs[i]
  
  # ------------------------------------------------------------
  # read in eQTL results
  # ------------------------------------------------------------
  
  pvs = read.table(inputi, header=TRUE, sep="\t", as.is=TRUE)
  dim(pvs)
  pvs[1:2,]
  
  summary(pvs$p.value)
  
  pvs = pvs[pvs$p.value < 1e-50,]
  dim(pvs)
  gc()
  
  summary(pvs$p.value)
  
  table(pvs$p.value < 1e-60)
  table(pvs$p.value < 1e-70)
  table(pvs$p.value < 1e-80)
  table(pvs$p.value < 1e-90)
  table(pvs$p.value < 1e-100)
    
  # ------------------------------------------------------------
  # plot it after removing those possible problematic probes
  # ------------------------------------------------------------
  
  table(pvs$SNP  %in% infMd$Name)
  pvs = pvs[which(pvs$SNP  %in% infMd$Name),]
  
  geneID   = match(pvs$gene, infoE$gene)
  markerID = match(pvs$SNP,  infoM$Composite.Element.REF)
  scores   = pvs$p.value
  scuts    = c(1e-60, 1e-70, 1e-80, 1e-100)
  cols     = c("green", "blue", "red", "black")
  eChr     = gsub("chr", "", infoE$chr)
  ePos     = 0.5*(infoE$start + infoE$end)
  mChr     = infoM$Chromosome
  mPos     = infoM$Genomic_Coordinate
  chroms   = 1:22
  
  outputi1 = sub(".png", "_rm_bad_m.png", outputi, fixed=TRUE)
  png(sprintf("../figures2/%s", outputi1), width=7.5, height=9,
  res=200, units="in")
  
  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab="Methylation Probe Location",
  ylab="Transcript Location", plot.hotspots=TRUE,
  hotspots.cut=10, score.type="p-value")
  
  dev.off()
  
}

q(save="no")
