
library(data.table)

# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# ------------------------------------------------------------------------
# check those associations between gene expression and DNA methylation
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# read in gene expression and methylation data
# ------------------------------------------------------------

datE = fread("expression_log_TReC_rm_cn_qnorm.txt", sep="\t", 
                   header=TRUE)
dim(datE)
datE[1:2,1:5]

datM = fread("methylation_mvalue_qnorm.txt", sep="\t", 
             header=TRUE)
dim(datM)
datM[1:2,1:5]

# ------------------------------------------------------------
# read in gene expression and methylation information
# ------------------------------------------------------------

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

infoE$chr = gsub("chr", "", infoE$chr)
dim(infoE)
infoE[1:2,]
table(infoE$chr, useNA="ifany")

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

pvs = read.table("expression_vs_methylation_rm_cn_pam50_qnorm.txt",
                  header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

summary(pvs$p.value)
table(pvs$SNP %in% infMd$Name)

pvs = pvs[which(pvs$SNP %in% infMd$Name),]
dim(pvs)
pvs[1:2,]
summary(pvs$p.value)

table(order(pvs$p.value) == 1:nrow(pvs))

table(pvs$p.value < 1e-50)
table(pvs$p.value < 1e-60)

# ------------------------------------------------------------------------
# check local EM pairs
# ------------------------------------------------------------------------

iM = infoM[match(pvs$SNP,  infoM$Composite.Element.REF),]
iE = infoE[match(pvs$gene, infoE$gene),]

dim(iM)
dim(iE)

iM[1:5,]
iE[1:5,]

table(iM$Chromosome == iE$chr, useNA="ifany")
table(iE$start < iE$end)

iE$pos = iE$start
iE$pos[which(iE$strand == "-")] = iE$end[which(iE$strand == "-")]

wLocal = which(iM$Chromosome == iE$chr & abs(iE$pos - iM$Genomic_Coordinate < 1e6))
length(wLocal)

pvs = pvs[wLocal,]
iM  = iM[wLocal,]
iE  = iE[wLocal,]

dim(pvs)
pvs[1:2,]
summary(pvs$p.value)

# ------------------------------------------------------------------------
# read in eQTL results after controlling for purities
# ------------------------------------------------------------------------

pvs7 = read.table("expression_vs_methylation_rm_cn_PCs7_qnorm.txt",
header=TRUE, sep="\t", as.is=TRUE)
dim(pvs7)
pvs7[1:2,]

summary(pvs7$p.value)

table(pvs7$SNP %in% infMd$Name)

pvs7 = pvs7[which(pvs7$SNP %in% infMd$Name),]
dim(pvs7)
pvs7[1:2,]
summary(pvs7$p.value)

table(order(pvs7$p.value) == 1:nrow(pvs7))

table(pvs7$p.value < 1e-20)
table(pvs7$p.value < 1e-30)
table(pvs7$p.value < 1e-50)

key1 = paste(pvs$gene, pvs$SNP, sep=":")
key7 = paste(pvs7$gene, pvs7$SNP, sep=":")

length(key7)
table(key1 %in% key7)

# ------------------------------------------------------------------------
# check local EM pairs after controlling for purities
# ------------------------------------------------------------------------

iM7 = infoM[match(pvs7$SNP,  infoM$Composite.Element.REF),]
iE7 = infoE[match(pvs7$gene, infoE$gene),]

dim(iM7)
dim(iE7)

iM7[1:5,]
iE7[1:5,]

table(iM7$Chromosome == iE7$chr, useNA="ifany")

iE7$pos = iE7$start
iE7$pos[which(iE7$strand == "-")] = iE7$end[which(iE7$strand == "-")]

wLocal7 = which(iM7$Chromosome == iE7$chr & abs(iE7$pos - iM7$Genomic_Coordinate < 1e6))
length(wLocal7)

pvs7 = pvs7[wLocal7,]
iM7  = iM7[wLocal7,]
iE7  = iE7[wLocal7,]

dim(pvs7)
pvs7[1:2,]
summary(pvs7$p.value)

table(pvs7$p.value < 1e-40)

# ------------------------------------------------------------------------
# check consistency between results contorlling purity or not
# ------------------------------------------------------------------------

key7 = paste(pvs7$gene, pvs7$SNP, sep=":")
table(key1 %in% key7)
table(key7 %in% key1)

keys  = union(key1, key7)
length(keys)
pval1 = rep(1e-36, length(keys))
pval7 = rep(1e-8,  length(keys))

pval1[match(key1, keys)] = pvs$p.value
pval7[match(key7, keys)] = pvs7$p.value

png("../figures2/EM_relation_with_purity_or_not_qnorm.png", width=7, height=3.5,
    res=400, units="in")

par(mfrow=c(1,2), mar=c(5.7,4,1.2,1), bty="n", cex=0.9)
plot(log10(-log10(pval1)), log10(-log10(pval7)), pch=20,
    col=rgb(1, 0.1, 0.1, 0.5),  xaxt="n", yaxt="n", cex=0.8,
    xlab="-log10(pval w/o purity correction)",
    ylab="-log10(pval with purity correction)")

abline(0, 1, lwd=2, col="darkgrey")
plbs = c(40, 60, 100, 200)
axis(side=1, at=log10(plbs), lab=plbs)

plbs = c(10, 20, 40, 100, 200)
axis(side=2, at=log10(plbs), lab=plbs)

# abline(h=log10(30), col="darkgrey", lty=5, lwd=1.5)

table(pval7 < 1e-40)
table(pval7 < 1e-30)

table(pval1 <= 1e-40, pval7 <= 1e-40)
table(pval1 <= 1e-40, pval7 <= 1e-30)

pval71 = -log10(pval7[which(pval1 <= 1e-40)])
length(pval71)
dim(pvs)
summary(pval71)

h1 = hist(pval71, breaks=c(5, 10, 20, 30, 40, 300), plot=FALSE)
mp = barplot(rev(h1$counts), ylim=c(0, max(h1$counts)*1.2))

txt = c("(0, 1e-40]", "(1e-40,1e-30]", "(1e-30,1e-20]", "(1e-20,1e-10]", ">1e-10")
text(mp, par("usr")[3], labels = txt, srt = 45, adj = c(1.1,1.1),
  xpd = TRUE, cex=.9)
text(mp, rev(h1$counts), labels=rev(h1$counts), pos=3)
mtext("pval w/o purity correction < 1e-40")
mtext("pval with purity correction", side=1, line=4.5)

dev.off()

# ------------------------------------------------------------------------
# check a few examples, 3 pairs with small p-value w/o purity correction
# ------------------------------------------------------------------------

cvx = read.table("cov_EM_with_pam50.txt", header=TRUE, sep="\t")
dim(cvx)
cvx[1:2,1:5]
X = t(data.matrix(cvx[,-1]))
dim(X)
X[1:2,1:5]

cvx = read.table("cov_EM_with_PCs7_qnorm.txt", header=TRUE, sep="\t")
dim(cvx)
cvx[1:2,1:5]
X2 = t(data.matrix(cvx[,-1]))
dim(X2)
X2[1:2,1:5]

table(rownames(X) == names(datE)[-1])

set.seed(2007)

table(pval1 <= 1e-80, pval7 >= 1e-10)
ww3 = sample(which(pval1 <= 1e-80 & pval7 >= 1e-10), 3)
keys[ww3]
ks = strsplit(keys[ww3], split=":")

pdf("../figures2/examples_local_EM_with_wo_purity_correction.pdf", 
    width=7.5, height=7.5)

par(mfrow=c(3,3), mar=c(5,4,1,1), bty="n", cex=0.9)

for(k in 1:3){
  gene1 = ks[[k]][1]
  cpg1  = ks[[k]][2]
  
  md1   = as.numeric(datM[which(datM$id == cpg1),-1])
  ed1   = as.numeric(datE[which(datE$id == gene1),-1])
  
  plot(md1, ed1, xlab=cpg1, ylab=gene1, cex=0.2)
  
  md2 = ed2 = rep(NA, length(md1))
  nNA = which((! is.na(md1)) & (! is.na(ed1)))
  
  md2[nNA] = lm(md1[nNA] ~ X[nNA,])$resid
  ed2[nNA] = lm(ed1[nNA] ~ X[nNA,])$resid
  plot(md2, ed2, xlab=cpg1, ylab=gene1, cex=0.2)
  
  md3 = ed3 = rep(NA, length(md1))
  nNA = which((! is.na(md1)) & (! is.na(ed1)))
  
  md3[nNA] = lm(md1[nNA] ~ X2[nNA,])$resid
  ed3[nNA] = lm(ed1[nNA] ~ X2[nNA,])$resid
  plot(md3, ed3, xlab=cpg1, ylab=gene1, cex=0.2)
}

dev.off()

q(save = "no")

