
# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# ------------------------------------------------------------
# read in tumor purity data
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

ff0 = "patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

table(is.na(emInfo$abs_purity))

# ------------------------------------------------------------------------------
# Read in eQTL results after removing copy number effects.
# ------------------------------------------------------------------------------

pvsCN = read.table("expression_vs_methylation_rm_cn_pam50_qnorm.txt",
                   header = TRUE, sep = "\t", as.is = TRUE)
dim(pvsCN)
pvsCN[1:2, ]

summary(pvsCN$p.value)

table(pvsCN$SNP  %in% infMd$Name)
pvsCN = pvsCN[which(pvsCN$SNP  %in% infMd$Name),]
dim(pvsCN)
pvsCN[1:2, ]
summary(pvsCN$p.value)

table(pvsCN$p.value < 1e-50)
table(pvsCN$p.value < 1e-60)
table(pvsCN$p.value < 1e-65)
table(pvsCN$p.value < 1e-70)

# Want to keep ~100-200k (methylation probe, gene expression) pairs.
pvs2kp = pvsCN[pvsCN$p.value < 1e-60, ]
dim(pvs2kp)

# ------------------------------------------------------------------------------
# Read in gene expression data and methylation information.
# ------------------------------------------------------------------------------

infoE = read.table("expression_info.txt",
                   sep = "\t", header = TRUE, as.is = TRUE)
dim(infoE)
infoE[1:2, ]

table(infoE$chr, useNA = "ifany")

infoM = read.table("methylation_info.txt",
                   sep = "\t", header = TRUE, as.is = TRUE)
dim(infoM)
infoM[1:2, ]
table(infoM$Chromosome, useNA = "ifany")

infoE$chr = gsub("chr", "", infoE$chr)
dim(infoE)
infoE[1:2, ]
table(infoE$chr, useNA = "ifany")

# ------------------------------------------------------------------------------
# Select pairs of (expression, methylation) so that for each pair
# expression and methylation are from different chromosomes.
# ------------------------------------------------------------------------------

iM = infoM[match(pvs2kp$SNP,  infoM$Composite.Element.REF),]
iE = infoE[match(pvs2kp$gene, infoE$gene),]

dim(iM)
dim(iE)

iM[1:5, ]
iE[1:5, ]

table(iM$Chromosome == iE$chr, useNA = "ifany")

pvs2kp = pvs2kp[which(iM$Chromosome != iE$chr),]
dim(pvs2kp)
pvs2kp[1:5, ]

length(unique(pvs2kp$gene))
length(unique(pvs2kp$SNP))

sort(table(pvs2kp$gene), decreasing = TRUE)[1:50]
sort(table(pvs2kp$SNP),  decreasing = TRUE)[1:50]

# ------------------------------------------------------------------------------
# Read in gene expression data and methylation data.
# ------------------------------------------------------------------------------

datE = read.table(file = "expression_log_TReC_qnorm.txt", sep = "\t",
                  header = TRUE, as.is = TRUE)
dim(datE)
datE[1:2, 1:5]

datM = read.table(file = "methylation_mvalue_qnorm.txt", sep = "\t",
                  header = TRUE, as.is = TRUE)
dim(datM)
datM[1:2, 1:5]

table(names(datE) == names(datM))

# datM is too big to do analysis below; restrict to only the relevant probes.
probes_to_keep = unique(pvsCN$SNP)

dim(pvsCN)
length(probes_to_keep)

dim(datM)
datM = datM[which(datM$id %in% probes_to_keep), ]
dim(datM)

# ------------------------------------------------------------------------------
# Take residuals.
# ------------------------------------------------------------------------------

pDatE = data.matrix(datE[, -1])
pDatM = data.matrix(datM[, -1])

cDat = read.table(file="cov_EM_with_pam50.txt", sep = "\t",
                  header = TRUE, as.is = TRUE)

dim(cDat)
cDat[1:5, 1:5]
cDat$id

X = data.matrix(t(cDat[, -1]))
dim(X)
X[1:5, 1:5]

X = cbind(rep(1, nrow(X)), X)
dim(X)
X[1:5, 1:5]

H = X %*% solve(t(X) %*% X) %*% t(X)

date()
rDatE = pDatE %*% (diag(nrow(H)) - H)
date()
dim(rDatE)

rDatM = pDatM %*% (diag(nrow(H)) - H)
date()
dim(rDatM)

# ------------------------------------------------------------------------------
# Prepare data of expression and methylation.
# ------------------------------------------------------------------------------

rDatE = scale(t(rDatE))
rDatM = scale(t(rDatM))

dim(rDatE)
dim(rDatM)

summary(apply(rDatE, 2, mean))
summary(apply(rDatE, 2, sd))

summary(apply(rDatM, 2, mean, na.rm = TRUE))
summary(apply(rDatM, 2, sd,   na.rm = TRUE))

rDatE1 = t(rDatE[, match(pvs2kp$gene, datE$id)])
rDatM1 = t(rDatM[, match(pvs2kp$SNP,  datM$id)])

dim(rDatE1)
dim(rDatM1)

rDatE1[1:5, 1:5]
rDatM1[1:5, 1:5]

# ------------------------------------------------------------------------------
# Plot the first 100 pairs and bottom.
# ------------------------------------------------------------------------------

pdf("../figures2/EM_1e-60_rm_cn_pcs_resid_ex_top100_qnorm.pdf",
    width = 5, height = 5)
par(mar = c(5,4,1,1), bty = "n")

# Some rows (probes) are problematic (all NA values). Skip such probes when
# printing out the top 100.
j = 0
for (i in 1:100) {
  check_next_row = TRUE
  while (check_next_row & j < min(nrow(rDatE1), nrow(rDatM1))) {
    j = j + 1
    ei = rDatE1[j, ]
    mi = rDatM1[j, ]
    if (length(which(is.na(ei))) == 0 &
        length(which(is.na(mi))) == 0) {
      check_next_row = FALSE
    } else {
      print(paste("Skipping row ", as.character(j), ": Num NA in ei: ",
                  as.character(length(which(is.na(ei)))), ", Num NA in mi: ",
                  as.character(length(which(is.na(mi)))), sep = ""))
    }
  }
  
  plot(mi, ei, xlab = pvs2kp$SNP[j], ylab = pvs2kp$gene[j])
}

dev.off()


pdf("../figures2/EM_1e-60_rm_cn_pcs_resid_ex_bottom100_qnorm.pdf",
    width = 5, height = 5)
par(mar = c(5, 4, 1, 1), bty = "n")

j = nrow(rDatE1) + 1
for (i in 1:100) {
  check_next_row = TRUE
  while (check_next_row & j > 0) {
    j = j - 1
    ei = rDatE1[j, ]
    mi = rDatM1[j, ]
    if (length(which(is.na(ei))) == 0 &
        length(which(is.na(mi))) == 0) {
      check_next_row = FALSE
    } else {
      print(paste("Skipping row ", as.character(j), ": Num NA in ei: ",
                  as.character(length(which(is.na(ei)))), ", Num NA in mi: ",
                  as.character(length(which(is.na(mi)))), sep = ""))
    }
  }
  
  plot(mi, ei, xlab = pvs2kp$SNP[j], ylab = pvs2kp$gene[j])
}

dev.off()

# ------------------------------------------------------------------------------
# Calculate PCs of each Pair.
# ------------------------------------------------------------------------------

signs = rep(1, nrow(pvs2kp))
signs[which(pvs2kp$beta < 0)] = -1
table(signs)

rDatM1[1:9, 1:5]
signs[1:9]
rDatM1 = rDatM1*signs
rDatM1[1:9, 1:5]

PCs = rDatE1 + rDatM1
PCs = PCs * sqrt(2) / 2

nna = rowSums(is.na(PCs))
table(nna)

PCs = PCs[which(nna == 0), ]
dim(PCs)

covdatR1 = t(PCs) %*% PCs / nrow(PCs)
dim(covdatR1)
prdatR1  = eigen(covdatR1)

prdatR1$values[1:20]
prdatR1$values[1:20]/sum(prdatR1$values)
cumsum(prdatR1$values)[1:20]/sum(prdatR1$values)

# ------------------------------------------------------------------------------
# Make the plot.
# ------------------------------------------------------------------------------

pdf("../figures2/ECM_pcs_qnorm.pdf",
   width = 4, height = 4)
par(mfrow = c(1, 1), mar = c(5, 4, 3, 1), bty = "n")
barplot(prdatR1$values[1:10], names.arg = 1:10)
dev.off()

pcs = prdatR1$vectors[, 1:50]
colnames(pcs) = paste("PC", 1:50, sep = "")

write.table(pcs, file = "ECM_pcs50_eigen_vectors_qnorm.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------
# make the plot again
# ------------------------------------------------------------------------

pdf("../figures2/purity_ECM_pam50_vs_purity_ABSOLUTE1_qnorm.pdf",
width=6, height=6)
par(mfrow=c(2,2), mar=c(5,4,3,1), bty="n", cex=0.9, cex.main=1, font.main=1)
barplot(prdatR1$values[1:10], names.arg=1:10, xlab="index of PCs",
ylab="Eigen-value")
for(j in 1:3){
  ct1 = cor.test(prdatR1$vectors[,j], emInfo$abs_purity)
  mm  = sprintf("corr=%.2f, p-value=%.1e", ct1$estimate, ct1$p.value)
  plot(prdatR1$vectors[,j], emInfo$abs_purity,
  xlab=sprintf("E-M PC %d", j),
  ylab="purity by ABSOLUTE", main=mm, cex=0.8)
}
dev.off()


pdf("../figures2/purity_ECM_pam50_vs_purity_ABSOLUTE2_qnorm.pdf",
width=6, height=6)
par(mfrow=c(2,2), mar=c(5,4,3,1), bty="n", cex=0.9, cex.main=1, font.main=1)
for(j in 4:7){
  ct1 = cor.test(prdatR1$vectors[,j], emInfo$abs_purity)
  mm  = sprintf("corr=%.2f, p-value=%.1e", ct1$estimate, ct1$p.value)
  plot(prdatR1$vectors[,j], emInfo$abs_purity,
  xlab=sprintf("E-M PC %d", j),
  ylab="purity by ABSOLUTE", main=mm, cex=0.8)
}
dev.off()

# ------------------------------------------------------------------------------
# Output new covariate data with one PC.
# ------------------------------------------------------------------------------

cDat = read.table(file="cov_EM_with_pam50.txt", sep = "\t",
                  header = TRUE, as.is = TRUE)

dim(cDat)
cDat[1:5, 1:5]
cDat$id

pt1 = data.frame(id = "pcs_ECM_PC1", t(pcs[, 1]))
names(pt1) = names(cDat)
dim(pt1)

cDat1 = rbind(cDat, pt1)
dim(cDat1)
cDat1[(nrow(cDat1) - 2):nrow(cDat1), 1:5]
cDat1$id

write.table(cDat1, file = "cov_EM_with_PCs1_qnorm.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# Output new covariate data with two PCs.
# ------------------------------------------------------------------------------

pt1 = data.frame(id = "pcs_ECM_PC2", t(pcs[, 2]))
names(pt1) = names(cDat)
dim(pt1)

cDat2 = rbind(cDat1, pt1)
dim(cDat2)
cDat2[(nrow(cDat2) - 2):nrow(cDat2), 1:5]
cDat2$id

write.table(cDat2, file = "cov_EM_with_PCs2_qnorm.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# Output new covariate data with three PCs.
# ------------------------------------------------------------------------------

pt1 = data.frame(id = "pcs_ECM_PC3", t(pcs[, 3]))
names(pt1) = names(cDat)
dim(pt1)

cDat3 = rbind(cDat2, pt1)
dim(cDat3)
cDat3[(nrow(cDat3) - 2):nrow(cDat3), 1:5]
cDat3$id

write.table(cDat3, file = "cov_EM_with_PCs3_qnorm.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# Output new covariate data with first 5 PCs.
# ------------------------------------------------------------------------------

pt1 = data.frame(id = "pcs_ECM_PC4", t(pcs[, 4]))

names(pt1) = names(cDat)
dim(pt1)

cDat4 = rbind(cDat3, pt1)

dim(cDat4)
cDat4[(nrow(cDat4) - 3):nrow(cDat4), 1:5]
cDat4$id

write.table(cDat4, file = "cov_EM_with_PCs4_qnorm.txt",
append = FALSE, quote = FALSE, sep = "\t",
row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------------------------
# Output new covariate data with first 5 PCs.
# ------------------------------------------------------------------------------

pt1 = data.frame(id = "pcs_ECM_PC5", t(pcs[, 5]))

names(pt1) = names(cDat)
dim(pt1)

cDat5 = rbind(cDat4, pt1)

dim(cDat5)
cDat5[(nrow(cDat5) - 3):nrow(cDat5), 1:5]
cDat5$id

write.table(cDat5, file = "cov_EM_with_PCs5_qnorm.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# Output new covariate data with first 6 PCs.
# ------------------------------------------------------------------------------

pt1 = data.frame(id = "pcs_ECM_PC6", t(pcs[, 6]))

names(pt1) = names(cDat)
dim(pt1)

cDat6 = rbind(cDat5, pt1)

dim(cDat6)
cDat6[(nrow(cDat6) - 3):nrow(cDat6), 1:5]
cDat6$id

write.table(cDat6, file = "cov_EM_with_PCs6_qnorm.txt",
append = FALSE, quote = FALSE, sep = "\t",
row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------------
# Output new covariate data with first 7 PCs.
# ------------------------------------------------------------------------------

pt1 = data.frame(id = "pcs_ECM_PC7", t(pcs[, 7]))

names(pt1) = names(cDat)
dim(pt1)

cDat7 = rbind(cDat6, pt1)

dim(cDat7)
cDat7[(nrow(cDat7) - 3):nrow(cDat7), 1:5]
cDat7$id

write.table(cDat7, file = "cov_EM_with_PCs7_qnorm.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

q(save="no")

