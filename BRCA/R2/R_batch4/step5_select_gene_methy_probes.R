
# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

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
# select top (expression, DNA methylation) pairs
# ------------------------------------------------------------------------

pvs = pvs[which(pvs$p.value < 1e-60),]
dim(pvs)

length(unique(pvs$SNP))
length(unique(pvs$gene))

summary(as.numeric(table(pvs$SNP)))
summary(as.numeric(table(pvs$gene)))

sort(table(pvs$SNP),  decreasing=TRUE)[1:20]
sort(table(pvs$gene), decreasing=TRUE)[1:20]

# ------------------------------------------------------------
# read in gene expression data and methylation information
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
# select pairs of (expression, methylation) so that for each pair
# expression and methylation are from different chromosomes
# ------------------------------------------------------------------------

iM = infoM[match(pvs$SNP,  infoM$Composite.Element.REF),]
iE = infoE[match(pvs$gene, infoE$gene),]

dim(iM)
dim(iE)

iM[1:5,]
iE[1:5,]

table(iM$Chromosome == iE$chr, useNA="ifany")

pvs2kp = pvs[which(iM$Chromosome != iE$chr),]
dim(pvs2kp)
pvs2kp[1:5,]

length(unique(pvs2kp$gene))
length(unique(pvs2kp$SNP))

# ------------------------------------------------------------------------
# obtain medain p-value and frequencies
# ------------------------------------------------------------------------

medPGe = tapply(pvs2kp$p.value, pvs2kp$gene, median)
tblGe  = table(pvs2kp$gene)

freqGe = as.numeric(tblGe)
medPGe = -log10(medPGe)

medPMe = tapply(pvs2kp$p.value, pvs2kp$SNP, median)
tblMe  = table(pvs2kp$SNP)

freqMe = as.numeric(tblMe)
medPMe = -log10(medPMe)

summary(freqGe)
summary(freqMe)

sort(tblGe, decreasing=TRUE)[1:20]
sort(tblMe, decreasing=TRUE)[1:20]

# ------------------------------------------------------------------------
# save results
# ------------------------------------------------------------------------

table(names(medPGe) == names(tblGe))
table(names(medPMe) == names(tblMe))

gene1e60 = data.frame(tblGe, medPGe)
meth1e60 = data.frame(tblMe, medPMe)

dim(gene1e60)
gene1e60[1:5,]
dim(meth1e60)
meth1e60[1:5,]

names(gene1e60) = c("gene", "freqIn1e60", "medianPvalIn1e60")
names(meth1e60) = c("methyProbe", "freqIn1e60", "medianPvalIn1e60")

gene1e60 = gene1e60[order(gene1e60$freqIn1e60, decreasing=TRUE),]
meth1e60 = meth1e60[order(meth1e60$freqIn1e60, decreasing=TRUE),]

dim(gene1e60)
gene1e60[1:5,]
dim(meth1e60)
meth1e60[1:5,]

write.table(gene1e60, file = "gene1e60_qnorm.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE,
col.names = TRUE)

write.table(meth1e60, file = "meth1e60_qnorm.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE,
col.names = TRUE)

# ------------------------------------------------------------
# read in gene expression data and methylation data
# ------------------------------------------------------------

datE = read.table(file = "expression_log_TReC_rm_cn_qnorm.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datE)
datE[1:2,1:5]

datM = read.table(file = "methylation_mvalue_qnorm.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datM)
datM[1:2,1:5]

table(names(datE) == names(datM))

# ------------------------------------------------------------------------
# read in EM PCs
# ------------------------------------------------------------------------

cDat = read.table("cov_EM_with_PCs7_qnorm.txt", sep="\t",
                  header=TRUE, as.is=TRUE)
dim(cDat)

dim(cDat)
cDat[1:2,1:5]

X = t(data.matrix(cDat[,-1]))
dim(X)

cDat$id

PCs = X[,43:49]
dim(PCs)

# ------------------------------------------------------------------------
# for each gene, check its association with top 7 PCs
# ------------------------------------------------------------------------

pDatE = data.matrix(datE[,-1])

Xdat0 = cbind(rep(1,nrow(X)), X[,1:42])
dim(Xdat0)
Xdat0[1:5,1:5]

table(colSums(is.na(Xdat0)))

H0 = Xdat0 %*% solve(t(Xdat0) %*% Xdat0) %*% t(Xdat0)

residE0 = pDatE %*% (diag(nrow(H0)) - H0)
dim(residE0)
residE0[1:2,1:5]

nnaE = rowSums(is.na(residE0))
table(nnaE)

# ------------------------------------------------------------------------
# the above approach will give some NAs, we re-calculate residuals
# ------------------------------------------------------------------------

wna = which(nnaE > 0)

for(i in 1:length(wna)){
  yi = pDatE[wna[i],]
  l1 = lm(yi ~ X[,1:42])
  residE0[wna[i], which(!is.na(yi))] = l1$resid
}

nnaE = rowSums(is.na(pDatE))
table(nnaE)

nnaE = rowSums(is.na(residE0))
table(nnaE)

# ------------------------------------------------------------------------
# for each methylation, check its association with top 7 PCs
# ------------------------------------------------------------------------

pDatM   = data.matrix(datM[,-1])
residM0 = pDatM %*% (diag(nrow(H0)) - H0)

dim(residM0)
residM0[1:2,1:5]

nnaM = rowSums(is.na(residM0))
table(nnaM)

# ------------------------------------------------------------------------
# the above approach will give some NAs, we re-calculate residuals
# ------------------------------------------------------------------------

wna  = which(nnaM > 0)

for(i in 1:length(wna)){
  if(i %% 1000 == 0){ cat(i, date(), "\n") }
  
  yi = pDatM[wna[i],]
  l1 = lm(yi ~ X[,1:42])
  residM0[wna[i], which(!is.na(yi))] = l1$resid
}

nnaM = rowSums(is.na(pDatM))
table(nnaM)

nnaM = rowSums(is.na(residM0))
table(nnaM)

# ------------------------------------------------------------------------
# calculate correlation between PCs and expression
# ------------------------------------------------------------------------

corE = cor(t(residE0), PCs, use = "pair")
dim(corE)
summary(corE)

colSums(abs(corE) > 0.5)
colSums(abs(corE) > 0.4)
colSums(abs(corE) > 0.3)
colSums(abs(corE) > 0.2)

df2 = ncol(residE0) - ncol(Xdat0) -1

Fs = (corE^2)/((1 - corE^2)/df2)
pE = pf(Fs, 1, df2, lower.tail=FALSE)

# ------------------------------------------------------------------------
# calculate correlation between PCs and methylation
# ------------------------------------------------------------------------

corM = cor(t(residM0), PCs, use = "pair")
dim(corM)
summary(corM)

colSums(abs(corM) > 0.5)
colSums(abs(corM) > 0.4)
colSums(abs(corM) > 0.3)
colSums(abs(corM) > 0.2)

df2 = ncol(residM0) - ncol(Xdat0) -1

Fs = (corM^2)/((1 - corM^2)/df2)
pM = pf(Fs, 1, df2, lower.tail=FALSE)

# ------------------------------------------------------------------------
# double check p-vlaues
# ------------------------------------------------------------------------

ptest = matrix(NA, nrow=1000, ncol=7)

for(i in 1:1000){
  yi = pDatM[i,]
  l1 = lm(yi ~ X)
  s1 = summary(l1)
  ptest[i,] = s1$coef[44:50,4]
}

round(cor(ptest, pM[1:1000,]),2)

# ------------------------------------------------------------------------
# plot it
# ------------------------------------------------------------------------

pdf("../figures2/cor_EM_PCs_qnorm.pdf", width=7, height=7)

for(i in 1:7){
  par(mfrow=c(2,2))

  hist(corE[,i], main=paste("PC", i, "vs. expression"),
      xlab="Correlation", xlim=c(-1,1), breaks=seq(-1,1,by=0.05))
      
  hist(pE[,i], main=paste("PC", i, "vs. expression"),
      xlab="p-value", xlim=c(0,1), breaks=seq(0,1,by=0.02))

  hist(corM[,i], main=paste("PC", i, "vs. methylation"),
      xlab="Correlation", xlim=c(-1,1), breaks=seq(-1,1,by=0.05))

  hist(pM[,i], main=paste("PC", i, "vs. methylation"),
      xlab="p-value", xlim=c(0,1), breaks=seq(0,1,by=0.02))

}

dev.off()

# ------------------------------------------------------------------------
# write out correlation matirx
# ------------------------------------------------------------------------

rownames(corE) = datE$id
colnames(corE) = paste("PC", 1:7, sep="")
dim(corE)
corE[1:2,]

write.table(corE, file = "corE_qnorm.txt", append = FALSE, quote = FALSE,
            sep = "\t", row.names = TRUE, col.names = TRUE)

rownames(corM) = datM$id
colnames(corM) = paste("PC", 1:7, sep="")
dim(corM)
corM[1:2,]

write.table(corM, file = "corM_qnorm.txt", append = FALSE, quote = FALSE,
            sep = "\t", row.names = TRUE, col.names = TRUE)

q(save = "no")

