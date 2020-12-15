
setwd("~/research/TCGA/_Sun_MethyE/BRCA")

ff0    = "_data2/patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

colSums(is.na(emInfo))

# ------------------------------------------------------------
# read in gene expression data
# ------------------------------------------------------------

datE = read.table(file = "RNAseqV2/rawCounts_tumor_v2.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datE)
datE[1:2,1:5]

table(emInfo$expression_barcode %in% gsub(".", "-", names(datE), fixed=TRUE))

cnms = strsplit(names(datE), split=".", fixed=TRUE)
table(sapply(cnms, length))
cnms = matrix(unlist(cnms), byrow=TRUE, ncol=7)

dim(cnms)
cnms[1:5,]
length(unique(cnms[,3]))

names(datE) = cnms[,3]

table(emInfo$patient_id %in% names(datE))

datE = datE[,match(emInfo$patient_id, names(datE))]
dim(datE)
datE[1:2,1:5]

# ------------------------------------------------------------
# read in gene location information
# ------------------------------------------------------------

ff2   = "_data2/gene_info.txt"
infoE = read.table(ff2, sep = "\t", header = TRUE, as.is=TRUE)
dim(infoE)
infoE[1:2,]
length(unique(infoE$gene))

table(rownames(datE) %in% infoE$gene)

features = intersect(rownames(datE), infoE$gene)
table(infoE$gene  %in% features)
table(rownames(datE)  %in% features)

datE  = datE[match(features, rownames(datE)),]
infoE = infoE[match(features, infoE$gene),]

dim(datE)
dim(infoE)

infoE[1:5,]

table(rownames(datE) == infoE$gene)
table(infoE$chr)
table(infoE$strand)
table(is.na(infoE$ensembl))

# ------------------------------------------------------------
# read in gene expression data direclty obtained from bam files
# ------------------------------------------------------------

datBam = read.table(file = "~/research/TCGA/BRCA/data/gene_counts_EA.txt",
  sep = "\t", header=TRUE, as.is=TRUE)
dim(datBam)
datBam[1:2,1:5]

genes2use = intersect(rownames(datBam), infoE$ensembl)
length(genes2use)

datBam2 = datBam[match(genes2use, rownames(datBam)), ]
datE2   = datE[match(genes2use, infoE$ensembl),]

dim(datBam2)
dim(datE2)
samples = intersect(names(datBam2), names(datE2))
length(samples)

datBam3 = datBam2[, match(samples, names(datBam2))]
datE3   = datE2[, match(samples, names(datE2))]

dim(datBam3)
dim(datE3)

datBam3[1:5,1:5]
datE3[1:5,1:5]

crs1 = crs2 = rep(NA, ncol(datBam3))

pdf("figures2/compare_TReC.pdf", width=6, height=6)
par(mar=c(5,4,1,1), bty="n")

for(i in 1:ncol(datBam3)){
  crs1[i] = cor(datBam3[,i], datE3[,i])
  crs2[i] = cor(log(datBam3[,i]+1), log(datE3[,i]+1))
  plot(log10(datBam3[,i]+1), log10(datE3[,i]+1),
      xlab="log10 TReC our own counting",
      ylab="log10 TReC from data portal",
      main=names(datBam3)[i])
  abline(0, 1, lwd=2, col="grey")
}
dev.off()

summary(crs1)
summary(crs2)

# ------------------------------------------------------------
# find a cutoff to filter out low expressed genes
# ------------------------------------------------------------

datEA = data.matrix(datE)

rMin = apply(datEA, 1, min)
rMed = apply(datEA, 1, median)
r75  = apply(datEA, 1, quantile, probs=0.75)
r90  = apply(datEA, 1, quantile, probs=0.90)

summary(rMin)
summary(rMed)
summary(r75)
summary(r90)

cor(rMin, rMed)
cor(r75,  rMed)
cor(r90,  rMed)

pdf("figures2/expression_cts_summary.pdf", width=6, height=6)
par(mfrow=c(2,2), mar=c(5,4,1,1), bty="n")
hist(log10(1+rMin), xlab="log10(min + 1)", main="")
hist(log10(1+rMed), xlab="log10(median + 1)", main="")
hist(log10(1+r75),  xlab="log10(75 percentile + 1)", main="")
hist(log10(1+r90),  xlab="log10(90 percentile + 1)", main="")
dev.off()

summary(rMin[rMed >=10])
summary(rMed[r75 >=20])

table(rMed >=10)
table(r75 >=20)

w2kp = which(r75 >=20)

dim(datEA)
datEA = datEA[w2kp,]
dim(datEA)

dim(infoE)
infoE = infoE[w2kp,]
dim(infoE)

if(! all(rownames(datEA) == infoE$gene) ){
  stop("gene name mismatch\n")
}

# ------------------------------------------------------------
# normalize gene expression by read-depth
# ------------------------------------------------------------

tot = colSums(datEA)
s75 = apply(datEA, 2, quantile, prob=0.75)

cor(tot, s75)

pdf("figures2/expression_total_vs_75_percentile.pdf", width=4, height=4)
par(mar=c(5,4,1,1), bty="n")
plot(tot/1e6, s75/1000, xlab="total reads (million)",
ylab="75 percentile (thousand)", cex=0.5)
dev.off()

nDat = t(log10(t((datEA + 1))/s75))
dim(nDat)

# ------------------------------------------------------------
# Run PCA using gene expression data, check possible outlier
# these PCs do include many batch effect information
# ------------------------------------------------------------

datR14Pr = nDat - rowMeans(nDat, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

pdf("figures2/expression_PCs_log_TReC.pdf", width=6, height=6)
par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="", xlab="Index", ylab="Eigen-value")

par(mar=c(5,4,1,1))
subtypes = unique(na.omit(emInfo$pam50))
cols = c("red", "green", "blue", "purple", "orange")

legend("topright", legend=subtypes, col=cols, pch=1, bty="n")

plot(PC1, PC2,  bty="n", cex=0.8)
for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(PC1[wj], PC2[wj], col=cols[j], cex=0.8)
}

plot(PC1, PC3,  bty="n", cex=0.8)
for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(PC1[wj], PC3[wj], col=cols[j], cex=0.8)
}

plot(PC2, PC3,  bty="n", cex=0.8)
for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(PC2[wj], PC3[wj], col=cols[j], cex=0.8)
}

dev.off()

cor(prdatR1$vectors[,1:5], emInfo$abs_purity, use="pair")
summary(lm(PC1 ~ emInfo$pam50))
summary(lm(PC2 ~ emInfo$pam50))
summary(lm(PC3 ~ emInfo$pam50))

# ------------------------------------------------------------
# write out data and information
# ------------------------------------------------------------

dim(datEA)
datEA[1:2,1:5]

write.table(datEA, file = "_data2/expression_counts.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = TRUE,col.names = TRUE)

table(colnames(nDat) == emInfo$patient_id)
pDat = data.frame(id=rownames(nDat), nDat)
dim(pDat)
pDat[1:2,1:5]

write.table(pDat, file = "_data2/expression_log_TReC.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(infoE)
infoE[1:5,]

write.table(infoE, file = "_data2/expression_info.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# check the association between gene expression and batches
# ------------------------------------------------------------

samE = strsplit(emInfo$expression_barcode, split="-", fixed=TRUE)
table(sapply(samE, length))
samE = matrix(unlist(samE), byrow=TRUE, ncol=7)
dim(samE)
samE[1,]

table(samE[,1])
table(samE[,7])

samE = data.frame(samE[,2:6], stringsAsFactors=FALSE)
names(samE) = c("institution", "patientID", "type", "portion", "plate")
dim(samE)
samE[1:2,]

table(rowSums(is.na(nDat)))

pE    = matrix(NA, nrow=nrow(nDat), ncol=6)
rDatE = matrix(NA, nrow=nrow(nDat), ncol=ncol(nDat))
age   = emInfo$age_at_diagnosis
PC1   = emInfo$noHM_PC1
PC2   = emInfo$noHM_PC2
PC3   = emInfo$noHM_PC3

for(i in 1:nrow(nDat)){
  
  if(i %% 1000 == 0){
    cat(i, date(), "\n")
  }
  
  yi = nDat[i,]
  li = lm(yi ~ samE$institution + samE$plate + age + PC1 + PC2 + PC3)
  ai = anova(li)
  
  pE[i,] = ai[1:6,5]

  rDatE[i,] = li$residuals
}

pdf("figures2/expression_vs_batches_p-value.pdf", width=10.5, height=7)
par(mfrow=c(2,3), mar=c(5,4,2,1), bty="n")
hist(pE[,1], xlab="p-value", main="instituion")
hist(pE[,2], xlab="p-value", main="plate")
hist(pE[,3], xlab="p-value", main="age")

hist(pE[,4], xlab="p-value", main="PC1")
hist(pE[,5], xlab="p-value", main="PC2")
hist(pE[,6], xlab="p-value", main="PC3")

dev.off()

summary(pE)

# ------------------------------------------------------------
# Run PCA using gene expression data, check possible outlier
# these PCs do include many batch effect information
# ------------------------------------------------------------

datR14Pr = rDatE - rowMeans(rDatE, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

pdf("figures2/expression_PCs_log_TReC_rm_institute_plate_age_genoPC.pdf",
width=6, height=6)
par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="",
xlab="Index", ylab="Eigen-value")

par(mar=c(5,4,1,1))
subtypes = unique(na.omit(emInfo$pam50))
cols = c("red", "green", "blue", "purple", "orange")

legend("topright", legend=subtypes, col=cols, pch=1, bty="n")

plot(PC1, PC2,  bty="n", cex=0.8)
for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(PC1[wj], PC2[wj], col=cols[j], cex=0.8)
}

plot(PC1, PC3,  bty="n", cex=0.8)
for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(PC1[wj], PC3[wj], col=cols[j], cex=0.8)
}

plot(PC2, PC3,  bty="n", cex=0.8)
for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(PC2[wj], PC3[wj], col=cols[j], cex=0.8)
}

dev.off()

cor(prdatR1$vectors[,1:5], emInfo$abs_purity, use="pair")
summary(lm(PC1 ~ emInfo$pam50))
summary(lm(PC2 ~ emInfo$pam50))
summary(lm(PC3 ~ emInfo$pam50))

# ------------------------------------------------------------
# Run clustreing on expression data, check possible outlier
# ------------------------------------------------------------

cr1   = cor(rDatE)
summary(as.numeric(cr1))

dist1 = as.dist(1 - cr1)
h1 = hclust(dist1)

pdf("figures2/expression_hclust_log_TReC_rm_institute_plate_age_genoPC.pdf",
    width=8, height=4)
par(mar=c(5,4,1,1), las=1)
plot(h1, main="", labels=FALSE, cex=0.7)

for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50[h1$order] == sj)
  mtext("I", side=1, line=0, col=cols[j], at=wj, cex=1, font=2)
}

legend("topright", legend=subtypes, col=cols, pch="I", bty="n", cex=0.9)

dev.off()

# if the median correaltion of one sample with other samples is too low
# indicating this sample is an outlier and we shoudl remove it

Ds = apply(cr1, 1, median)
pdf("figures2/expression_D_statistics_rm_institute_plate_age_genoPC.pdf",
    width=8, height=4)
par(mar=c(5,4,1,1), mfrow=c(1,2), bty="n")
hist(Ds, xlab="D statistics", breaks=20, main="")

plot(Ds, emInfo$abs_purity, type="n", ylab="Purity")

for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(Ds[wj], emInfo$abs_purity[wj], col=cols[j], cex=0.8)
}

dev.off()



q(save = "no")
