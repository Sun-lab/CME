
library(gplots)

# ------------------------------------------------------------
# read in the information of the samples to be used
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA")

ff0    = "_data2/patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

# ------------------------------------------------------------
# read in copy number data
# ------------------------------------------------------------

datC = read.table(file = "CNV_SNP_Array/cn_data.txt", sep = "\t",
                  header = TRUE, as.is=TRUE)
dim(datC)
datC[1:2,1:5]

length(unique(names(datC)))
table(emInfo$patient_id %in% names(datC))

datC = datC[,match(emInfo$patient_id, names(datC))]
names(datC) = emInfo$patient_id

dim(datC)
datC[1:2,1:5]

# ------------------------------------------------------------
# remove probes with many NAs and summary mean/sds
# ------------------------------------------------------------

nna  = rowSums(is.na(datC))
summary(nna)
table(nna[nna > 0])

0.10*ncol(datC)
table(nna < 0.10*ncol(datC))

0.05*ncol(datC)
table(nna < 0.05*ncol(datC))

table(nna < 1)

w2kp = which(nna < 0.05*ncol(datC))
datC = datC[w2kp,]
dim(datC)
datC[1:2,1:5]

nna  = colSums(is.na(datC))
summary(nna)
table(nna[nna > 0])

avs = rowMeans(datC, na.rm=TRUE)
sds = apply(datC, 1, sd, na.rm=TRUE)

summary(sds)
summary(avs)

pdf("figures2/cn_values.pdf", width=7, height=3.5)

mat = matrix(c(1,2,3,3), ncol=2)
layout(mat)

par(mar=c(5, 4, 1, 1), bty="n")
hist(avs, xlab="mean values", main="")
hist(sds, xlab="sds", main="")
smoothScatter(avs, sds, xlab="mean values", ylab="sds")

dev.off()

# ------------------------------------------------------------
# read in gene location information
# ------------------------------------------------------------

ff0   = "_data2/gene_info.txt"
infoC = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE, comment.char="")
dim(infoC)
infoC[1:2,]
length(unique(infoC$gene))

table(rownames(datC) %in% infoC$gene)

infoC = infoC[match(rownames(datC), infoC$gene),]
dim(infoC)
infoC[1:5,]

table(rownames(datC) == infoC$gene)
table(infoC$chr, useNA="ifany")

# ------------------------------------------------------------
# Run PCA using copy number data, check possible outlier
# these PCs do include many batch effect information
# ------------------------------------------------------------

datC     = data.matrix(datC)
datR14Pr = datC - rowMeans(datC, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

pdf("figures2/cn_PCs.pdf", width=6, height=6)
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
cor(prdatR1$vectors[,1:5], emInfo$age_at_diagnosis, use="pair")

summary(lm(PC1 ~ emInfo$pam50))
summary(lm(PC2 ~ emInfo$pam50))
summary(lm(PC3 ~ emInfo$pam50))

# ------------------------------------------------------------
# order genes by locations
# ------------------------------------------------------------

table(infoC$chr, useNA="ifany")
numChr = gsub("chr", "", infoC$chr)
numChr[which(numChr == "X")] = "23"
numChr = as.numeric(numChr)
table(numChr, useNA="ifany")

od1   = order(numChr, infoC$start)
datC  = datC[od1,]
infoC = infoC[od1,]

# ------------------------------------------------------------
# write out data and information
# ------------------------------------------------------------

table(colnames(datC) == emInfo$patient_id)
table(rownames(datC) == infoC$gene)

pDat = data.frame(id=rownames(datC), datC)
dim(pDat)
pDat[1:2,1:5]

write.table(pDat, file = "_data2/cn_values.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(infoC)
infoC[1:5,]

write.table(infoC, file = "_data2/cn_info.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

################################################
##### Truncated for visualization purposes #####
################################################

summary(as.numeric(datC))

toplot = datC
toplot[which(datC > 1)]  = 1
toplot[which(datC < -1)] = -1

binnames   = infoC$chr
binlabel   = gsub("chr", "", binnames)
binlabel   = gsub("X", "23", binlabel)
difference = diff(as.numeric(binlabel))
seps       = which(difference != 0)+1
rowseps    = c(1, seps, length(binnames))
binlabel   = rep("", length(binnames))
middle     = round((c(1, seps) + c(seps, length(binnames)))/2)
binlabel[middle] = c(1:22, "X")


palette.gr.marray = colorRampPalette(c("blue", "white", "red"))
png("figures2/cn_heatmap.png", width = 10, height = 10, units="in", res=400)

h2 = heatmap.2(toplot, trace = "none", col = palette.gr.marray,
          Rowv = F, dendrogram = "column", key = T, symbreaks = T,
          sepwidth = c(0.02, 0.02), rowsep = rowseps, sepcolor = "gray30",
          labRow = binlabel, cexRow = 0.7, labCol = NA, cexCol = 1)
dev.off()

# ------------------------------------------------------------
# Run clustreing on cn data
# ------------------------------------------------------------

dist1 = dist(t(toplot))
h1    = hclust(dist1)

pdf("figures2/cn_hclust.pdf", width=8, height=4)
par(mar=c(5,4,1,1))
plot(h1, main="", labels=FALSE, cex=0.7, hang=-1)

for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50[h1$order] == sj)
  mtext("I", side=1, line=0, col=cols[j], at=wj, cex=1, font=2)
}

legend("topright", legend=subtypes, col=cols, pch="I", bty="n", cex=0.9)
dev.off()

# ------------------------------------------------------------
# obtain batch information for copy number data
# ------------------------------------------------------------

samC = strsplit(emInfo$cn_barcode, split="-", fixed=TRUE)
table(sapply(samC, length))

samC = matrix(unlist(samC), byrow=TRUE, ncol=7)
dim(samC)
samC[1:2,]

samC = samC[,2:6]
dim(samC)
samC[1:2,]

samC = data.frame(samC, stringsAsFactors=FALSE)
names(samC) = c("institution", "patientID", "type", "portion", "plate")
dim(samC)
samC[1:2,]

samM = strsplit(emInfo$methylation_file, split="-", fixed=TRUE)
table(sapply(samM, length))
samM = matrix(unlist(samM), byrow=TRUE, ncol=9)
dim(samM)
samM[1:2,]

samM = samM[,4:8]
dim(samM)
samM[1:2,]

samM = data.frame(samM, stringsAsFactors=FALSE)
names(samM) = c("institution", "patientID", "type", "portion", "plate")
dim(samM)
samM[1:2,]


table(samC$institution == samM$institution)
table(samC$patientID == samM$patientID)

sort(table(samC$plate))
sort(table(samM$plate))

length(table(samC$plate))
length(table(samM$plate))

## the plate informtion of copy number and methylation are almost the same
## so I will just use the plate informaiton from expression data
table(samC$plate, samM$plate)

# ------------------------------------------------------------
# check the association between copy number and batches
# ------------------------------------------------------------

table(rowSums(is.na(datC)))

pC    = matrix(NA, nrow=nrow(datC), ncol=6)
rDatC = matrix(NA, nrow=nrow(datC), ncol=ncol(datC))
age   = emInfo$age_at_diagnosis
PC1   = emInfo$noHM_PC1
PC2   = emInfo$noHM_PC2
PC3   = emInfo$noHM_PC3

for(i in 1:nrow(datC)){
  
  if(i %% 1000 == 0){
    cat(i, date(), "\n")
  }
  
  yi = datC[i,]
  li = lm(yi ~ samC$institution + samC$plate + age + PC1 + PC2 + PC3)
  ai = anova(li)
  
  pC[i,] = ai[1:6,5]
  
  rDatC[i,which(!is.na(yi))] = li$residuals
}

pdf("figures2/cn_vs_batches_p-value.pdf", width=10.5, height=7)
par(mfrow=c(2,3), mar=c(5,4,2,1), bty="n")
hist(pC[,1], xlab="p-value", main="instituion", breaks=20)
hist(pC[,2], xlab="p-value", main="plate", breaks=20)
hist(pC[,3], xlab="p-value", main="age", breaks=20)

hist(pC[,4], xlab="p-value", main="PC1", breaks=20)
hist(pC[,5], xlab="p-value", main="PC2", breaks=20)
hist(pC[,6], xlab="p-value", main="PC3", breaks=20)

dev.off()

summary(pC)

table(pC[,3] < 1e-4)
table(pC[,3] < 1e-5)
table(pC[,3] < 1e-6)

cbind(infoC, pC[,3:4])[which(pC[,3] < 1e-5),]

IL23R = datC[which(infoC$gene=="IL23R|149233"),]
summary(lm(IL23R ~ age))

pdf("figures2/cn_vs_age_IL23R.pdf", width=4, height=4)
par(mar=c(5,4,2,1), bty="n")
plot(age, IL23R, pch=20, col=rgb(1, 0.1, 0.1, 0.5))
dev.off()

# ------------------------------------------------------------
# Manhattan for age
# ------------------------------------------------------------

source("../shared_code/manhattan.R")

png("figures2/cn_vs_age.png", width=10, height=3.5, res=400, units="in")
par(mar=c(5,4,1,1), bty="n")
manhattan(pC[,3], gsub("chr", "", infoC$chr), infoC$start)
dev.off()

# ------------------------------------------------------------
# draw heatmap using residualized data
# ------------------------------------------------------------

summary(as.numeric(rDatC))

toplot = rDatC
toplot[which(rDatC > 1)]  = 1
toplot[which(rDatC < -1)] = -1

png("figures2/cn_heatmap_rm_institute_plate_age_genoPC.png",
    width = 10, height = 10, units="in", res=400)

heatmap.2(toplot, trace = "none", col = palette.gr.marray,
Rowv = F, dendrogram = "column", key = T, symbreaks = T,
sepwidth = c(0.02, 0.02), rowsep = rowseps, sepcolor = "gray30",
labRow = binlabel, cexRow = 0.7, labCol = NA, cexCol = 1)
dev.off()

# ------------------------------------------------------------
# Run clustreing on residualized cn data
# ------------------------------------------------------------

dist1 = dist(t(toplot))
h1    = hclust(dist1)

pdf("figures2/cn_hclust_rm_institute_plate_age_genoPC.pdf",
width=8, height=4)
par(mar=c(5,4,1,1))
plot(h1, main="", labels=FALSE, cex=0.7, hang=-1)

for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50[h1$order] == sj)
  mtext("I", side=1, line=0, col=cols[j], at=wj, cex=1, font=2)
}

legend("topright", legend=subtypes, col=cols, pch="I", bty="n", cex=0.9)
dev.off()

# ------------------------------------------------------------
# Run PCA using cn data, check possible outlier
# these PCs do include many batch effect information
# ------------------------------------------------------------

datR14Pr = rDatC - rowMeans(rDatC, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

pdf("figures2/cn_PCs_mvalue_rm_institute_plate_age_genoPCs.pdf",
width=6, height=6)
par(mar=c(5,4,1,1), mfrow=c(2,2))
barplot(prdatR1$values[1:20], main="",
xlab="Index", ylab="Eigen-value")

subtypes = unique(na.omit(emInfo$pam50))
cols = c("red", "green", "blue", "purple", "orange")

legend("topright", legend=subtypes, col=cols, pch=1, bty="n")

par(mar=c(5,4,1,1))

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
# if the median correaltion of one sample with other samples is too low
# indicating this sample is an outlier and we shoudl remove it
# ------------------------------------------------------------

cr1 = cor(rDatC, use="pair")
Ds  = apply(cr1, 1, median)
pdf("figures2/cn_D_statistics_rm_institute_plate_age_genoPC.pdf",
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

emInfo$Ds = Ds
emInfo[emInfo$Ds < quantile(Ds, 0.05), c(1:2,20:22)]


q(save = "no")
