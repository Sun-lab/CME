
require(gplots)

# ------------------------------------------------------------
# read in the information of the samples to be used
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA")

ff0    = "_data2/patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

# ------------------------------------------------------------
# read in DNA methylation data
# ------------------------------------------------------------

datM = read.table(file = "DNA_Methylation/mValues_tumor_v2.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datM)
datM[1:2,1:5]

ss1 = gsub(".", "-", names(datM), fixed=TRUE)
ss1 = gsub("-TCGA", "/TCGA", ss1, fixed=TRUE)

table(emInfo$methylation_barcode %in% ss1)

datM = datM[,match(emInfo$methylation_barcode, ss1)]
names(datM) = emInfo$patient_id

dim(datM)
datM[1:2,1:5]

# ------------------------------------------------------------
# remove probes with many NAs and summary mean/sds
# ------------------------------------------------------------

nna  = rowSums(is.na(datM))
summary(nna)
table(nna[nna > 0])

0.10*ncol(datM)
table(nna < 0.10*ncol(datM))

0.05*ncol(datM)
table(nna < 0.05*ncol(datM))

table(nna < 1)

w2kp = which(nna < 0.05*ncol(datM))
datM = datM[w2kp,]
dim(datM)
datM[1:2,1:5]

avs = rowMeans(datM, na.rm=TRUE)
sds = apply(datM, 1, sd, na.rm=TRUE)

summary(sds)
summary(avs)

pdf("figures2/methylation_mvalues.pdf", width=7, height=3.5)

mat = matrix(c(1,2,3,3), ncol=2)
layout(mat)

par(mar=c(5, 4, 1, 1), bty="n")
hist(avs, xlab="mean values", main="")
hist(sds, xlab="sds", main="")
smoothScatter(avs, sds, xlab="mean values", ylab="sds")

dev.off()

# ------------------------------------------------------------
# read in DNA methylation information
# ------------------------------------------------------------

ff0   = "DNA_Methylation/HM450_info.txt"
infoM = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE, comment.char="")
dim(infoM)
infoM[1:2,]
length(unique(infoM$Composite.Element.REF))

table(rownames(datM) %in% infoM$Composite.Element.REF)

infoM = infoM[match(rownames(datM), infoM$Composite.Element.REF),]
dim(infoM)
infoM[1:5,]

table(rownames(datM) == infoM$Composite.Element.REF)
table(infoM$Chromosome, useNA="ifany")

infoM[is.na(infoM$Chromosome),]

numChr = infoM$Chromosome
numChr[infoM$Chromosome=="X"] = "23"
numChr[infoM$Chromosome=="Y"] = "24"
numChr = as.numeric(numChr)
table(numChr)

od1   = order(numChr, infoM$Genomic_Coordinate)
datM  = datM[od1,]
infoM = infoM[od1,]

# ------------------------------------------------------------
# plot heatmap of DNA methylation
# ------------------------------------------------------------

toplot = data.matrix(exp(datM)/(1 + exp(datM)))
sds    = apply(toplot, 1, sd, na.rm=TRUE)
summary(sds)
table(sds > 0.2)
table(sds > 0.15)

w2kp   = which(sds > 0.2 & !is.na(infoM$Chromosome) & infoM$Chromosome!="Y")
toplot = toplot[w2kp,]
infoC  = infoM[w2kp,]

dim(toplot)
toplot[1:2,1:5]

dim(infoC)
infoC[1:2,]

table(rownames(toplot) == infoC$Composite.Element.REF)
summary(as.numeric(toplot))

png("figures2/methylation_heatmap.png", width = 10, height = 10, units="in", res=400)

palette.gr.marray = colorRampPalette(c("blue", "white", "red"))
h2 = heatmap.2(toplot, trace = "none", col = palette.gr.marray,
              Rowv = T, dendrogram = "both", key = T, symbreaks = F,
              symkey = F, labRow = NA, labCol = NA, cexCol = 1)

dev.off()

h2$carpet = NULL
save(h2, file = "_data2/methylation_cluster.RData")

# ------------------------------------------------------------
# plot density of DNA methylation
# ------------------------------------------------------------

hc1 = as.hclust(h2$colDendrogram)
c1  = cutree(hc1, k=6)
table(names(c1) == colnames(toplot))

png("figures2/methylation_density_by_cluster.png", width = 12,
    height = 8, units="in", res=400)
par(mfrow=c(2,3), mar=c(5, 4, 1, 1), bty="n")
cols = rainbow(6)
for(k in 1:6){
  plot(c(0,1), c(0,5), type="n", xlab="beta-value", ylab="density", main="")
  wk = which(c1==k)
  for(wkj in wk){
    lines(density(na.omit(toplot[,wkj])), col=cols[k])
  }
}

dev.off()

# ------------------------------------------------------------
# Run PCA using methylation data, check possible outlier
# these PCs do include many batch effect information
# ------------------------------------------------------------

datM     = data.matrix(datM)
datR14Pr = datM - rowMeans(datM, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

pdf("figures2/methylation_PCs_mvalue.pdf", width=6, height=6)
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

cor(prdatR1$vectors[,1:5], emInfo$absolute_extract_purity, use="pair")
cor(prdatR1$vectors[,1:5], emInfo$age_at_diagnosis, use="pair")

summary(lm(PC1 ~ emInfo$pam50))
summary(lm(PC2 ~ emInfo$pam50))
summary(lm(PC3 ~ emInfo$pam50))

# ------------------------------------------------------------
# write out data and information
# ------------------------------------------------------------

table(colnames(datM) == emInfo$patient_id)

pDat = data.frame(id=rownames(datM), datM)
dim(pDat)
pDat[1:2,1:5]

write.table(pDat, file = "_data2/methylation_mvalue.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(infoM)
infoM[1:5,]

write.table(infoM, file = "_data2/methylation_info.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# check the association between DNA methylation and batches
# ------------------------------------------------------------

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

table(rowSums(is.na(datM)))

pM    = matrix(NA, nrow=nrow(datM), ncol=6)
rDatM = matrix(NA, nrow=nrow(datM), ncol=ncol(datM))
age   = emInfo$age_at_diagnosis
PC1   = emInfo$noHM_PC1
PC2   = emInfo$noHM_PC2
PC3   = emInfo$noHM_PC3

for(i in 1:nrow(datM)){
  
  if(i %% 1000 == 0){
    cat(i, date(), "\n")
  }
  
  yi = datM[i,]
  li = lm(yi ~ samM$institution + samM$plate + age + PC1 + PC2 + PC3)
  ai = anova(li)
  
  pM[i,] = ai[1:6,5]

  rDatM[i,which(!is.na(yi))] = li$residuals
}

pdf("figures2/methylation_vs_batches_p-value.pdf", width=10.5, height=7)
par(mfrow=c(2,3), mar=c(5,4,2,1), bty="n")
hist(pM[,1], xlab="p-value", main="instituion")
hist(pM[,2], xlab="p-value", main="plate")
hist(pM[,3], xlab="p-value", main="age")

hist(pM[,4], xlab="p-value", main="PC1")
hist(pM[,5], xlab="p-value", main="PC2")
hist(pM[,6], xlab="p-value", main="PC3")

dev.off()

summary(pM)

# ------------------------------------------------------------
# Run PCA using methylation data, check possible outlier
# these PCs do include many batch effect information
# ------------------------------------------------------------

datR14Pr = rDatM - rowMeans(rDatM, na.rm=TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]

PC1 =  prdatR1$vectors[,1]
PC2 =  prdatR1$vectors[,2]
PC3 =  prdatR1$vectors[,3]

pdf("figures2/methylation_PCs_mvalue_rm_institute_plate_age_genoPC.pdf",
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

cor(prdatR1$vectors[,1:5], emInfo$absolute_extract_purity, use="pair")

summary(lm(PC1 ~ emInfo$pam50))
summary(lm(PC2 ~ emInfo$pam50))
summary(lm(PC3 ~ emInfo$pam50))

# ------------------------------------------------------------
# Run clustreing on methylation data, check possible outlier
# ------------------------------------------------------------

cr1   = cor(rDatM, use="pair")
summary(as.numeric(cr1))

dist1 = as.dist(1 - cr1)
h1 = hclust(dist1)

pdf("figures2/methylation_hclust_mvalue_rm_institute_plate_age_genoPC.pdf",
    width=8, height=4)
par(mar=c(5,4,1,1))
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
pdf("figures2/methylation_D_statistics_rm_institute_plate_age_genoPC.pdf",
    width=8, height=4)
par(mar=c(5,4,1,1), mfrow=c(1,2), bty="n")
hist(Ds, xlab="D statistics", breaks=20, main="")
plot(Ds, emInfo$absolute_extract_purity, type="n", ylab="Purity")

for(j in 1:5){
  sj = subtypes[j]
  wj = which(emInfo$pam50 == sj)
  points(Ds[wj], emInfo$absolute_extract_purity[wj], col=cols[j], cex=0.8)
}
dev.off()

emInfo$Ds = Ds
emInfo[emInfo$Ds < quantile(Ds, 0.05), c(1:2,20:22)]


q(save = "no")
