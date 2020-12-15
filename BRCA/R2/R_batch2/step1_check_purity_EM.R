
# ------------------------------------------------------------
# read in tumor purity data
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

ff0 = "patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

table(is.na(emInfo$abs_purity))

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

pvs = read.table("expression_vs_methylation.txt",
                  header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

table(pvs$p.value < 1e-50)
table(pvs$p.value < 1e-60)
table(pvs$p.value < 1e-80)

pvs2kp = pvs[pvs$p.value < 1e-60,]
dim(pvs2kp)

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

iM = infoM[match(pvs2kp$SNP,  infoM$Composite.Element.REF),]
iE = infoE[match(pvs2kp$gene, infoE$gene),]

dim(iM)
dim(iE)

iM[1:5,]
iE[1:5,]

table(iM$Chromosome == iE$chr, useNA="ifany")

pvs2kp = pvs2kp[which(iM$Chromosome != iE$chr),]
dim(pvs2kp)
pvs2kp[1:5,]

length(unique(pvs2kp$gene))
length(unique(pvs2kp$SNP))

sort(table(pvs2kp$gene), decreasing=TRUE)[1:50]
sort(table(pvs2kp$SNP),  decreasing=TRUE)[1:50]

# ------------------------------------------------------------
# read in gene expression data and methylation data
# ------------------------------------------------------------

datE = read.table(file = "expression_log_TReC.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datE)
datE[1:2,1:5]

datM = read.table(file = "methylation_mvalue.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datM)
datM[1:2,1:5]

table(names(datE) == names(datM))

pDatE = data.matrix(datE[,-1])
pDatM = data.matrix(datM[,-1])

# ------------------------------------------------------------------------
# prepare data of expression and methylation
# ------------------------------------------------------------------------

pDatE = scale(t(pDatE))
pDatM = scale(t(pDatM))

dim(pDatE)
dim(pDatM)

summary(apply(pDatE, 2, mean))
summary(apply(pDatE, 2, sd))

summary(apply(pDatM, 2, mean, na.rm=TRUE))
summary(apply(pDatM, 2, sd,   na.rm=TRUE))

pDatE1 = t(pDatE[,match(pvs2kp$gene, datE$id)])
pDatM1 = t(pDatM[,match(pvs2kp$SNP,  datM$id)])

dim(pDatE1)
dim(pDatM1)

pDatE1[1:5,1:5]
pDatM1[1:5,1:5]

# ------------------------------------------------------------------------
# plot the first 100 pairs and bottom
# ------------------------------------------------------------------------

pdf("../figures2/EM_1e-30_ex_top100.pdf", width=5, height=5)
par(mar=c(5,4,1,1), bty="n")

for(i in 1:100){
  
  ei = pDatE1[i,]
  mi = pDatM1[i,]
  
  plot(mi, ei, xlab=pvs$SNP[i], ylab=pvs$gene[i])
}

dev.off()


pdf("../figures2/EM_1e-30_ex_bottom100.pdf", width=5, height=5)
par(mar=c(5,4,1,1), bty="n")

for(i in nrow(pDatE1):(nrow(pDatE1)-99)){
  
  ei = pDatE1[i,]
  mi = pDatM1[i,]
  
  plot(mi, ei, xlab=pvs$SNP[i], ylab=pvs$gene[i])
}

dev.off()

# ------------------------------------------------------------------------
# calculate PCs of each Pair
# ------------------------------------------------------------------------

signs = rep(1, nrow(pvs2kp))
signs[which(pvs2kp$beta < 0)] = -1
table(signs)

pDatM1[1:9,1:5]
signs[1:9]
pDatM1 = pDatM1*signs
pDatM1[1:9,1:5]

PCs = pDatE1 + pDatM1
PCs = PCs*sqrt(2)/2

nna = rowSums(is.na(PCs))
table(nna)

PCs = PCs[which(nna == 0),]
dim(PCs)

covdatR1 = t(PCs) %*% PCs / nrow(PCs)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]
cr1 = cor(prdatR1$vectors[,1:10], emInfo$abs_purity, use="pair")
cr1 = round(cr1, 4)

# ------------------------------------------------------------------------
# Some of these PCs are actually correlated with batch effects
# ------------------------------------------------------------------------

barcode = strsplit(emInfo$methylation_file, split="-")
table(sapply(barcode, length))

barcode = matrix(unlist(barcode), ncol=9, byrow=TRUE)

dim(barcode)
barcode[1:2,]

plate = barcode[,8]
table(plate)

site  = emInfo$tissue_source_site
sort(table(site))

stage = emInfo$ajcc_pathologic_tumor_stage
pam50 = emInfo$pam50

age = emInfo$age_at_diagnosis
table(is.na(age))

table(emInfo$Path.Score, useNA="ifany")

pvs1 = matrix(NA, nrow=10, ncol=5)

for(i in 1:10){
  yi = prdatR1$vectors[,i]
  li = lm(yi ~ site + plate + stage + age + pam50)
  ai = anova(li)
  
  pvs1[i,] = ai[1:5,5]
}

pvs1 = signif(pvs1,2)
colnames(pvs1) = c("site", "plate", "stage", "age", "pam50")

cbind(cr1, pvs1)

# ------------------------------------------------------------------------
# make the plot
# ------------------------------------------------------------------------

pdf("../figures2/purity_EM_vs_purity_ABSOLUTE.pdf",
    width=7, height=10.5)

par(mfrow=c(3,2), mar=c(5,4,3,1), bty="n")

barplot(prdatR1$values[1:10], names.arg=1:10)

for(j in 1:5){
  ct1 = cor.test(prdatR1$vectors[,j], emInfo$abs_purity)
  mm  = sprintf("corr=%.2f, p-value=%.1e", ct1$estimate, ct1$p.value)
  plot(prdatR1$vectors[,j], emInfo$abs_purity,
  xlab=sprintf("purity by EM, PC %d", j),
  ylab="purity by ABSOLUTE", main=mm)
}

dev.off()

pcs = prdatR1$vectors[,1:50]
colnames(pcs) = paste("PC", 1:50, sep="")

write.table(pcs, file = "EM_eigen_vectors.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)


q(save = "no")

