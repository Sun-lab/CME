
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

infoM = read.table("methylation_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoM)
infoM[1:2,]

table(infoM$Chromosome, useNA="ifany")

infoC = read.table("cn_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoC)
infoC[1:2,]
table(infoC$chr, useNA="ifany")

infoC$pos = 0.5*(infoC$start + infoC$end)
names(infoM)[3:4] = c("chr", "pos")
infoM$chr = paste("chr", infoM$chr, sep="")
dim(infoM)
infoM[1:2,]

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

pvs0 = read.table("methylation_vs_cn.txt", header=TRUE, sep="\t", as.is=TRUE)
dim(pvs0)
pvs0[1:2,]

summary(pvs0$p.value)
summary(pvs0$FDR)

table(pvs0$p.value < 1e-30)

pvs0$CN.chr = infoC$chr[match(pvs0$SNP, infoC$gene)]
pvs0$CN.pos = infoC$pos[match(pvs0$SNP, infoC$gene)]

pvs0$gene.chr = infoM$chr[match(pvs0$gene, infoM$Composite.Element.REF)]
pvs0$gene.pos = infoM$pos[match(pvs0$gene, infoM$Composite.Element.REF)]

table(pvs0$CN.chr == pvs0$gene.chr)
table(pvs0$CN.chr == pvs0$gene.chr)/nrow(pvs0)

# ------------------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------------------

pvs = read.table("methylation_vs_cn_with_ab_purity_pam50.txt",
                  header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

summary(pvs$p.value)
summary(pvs$FDR)

table(pvs$p.value < 1e-30)

pvs$CN.chr = infoC$chr[match(pvs$SNP, infoC$gene)]
pvs$CN.pos = infoC$pos[match(pvs$SNP, infoC$gene)]

pvs$gene.chr = infoM$chr[match(pvs$gene, infoM$Composite.Element.REF)]
pvs$gene.pos = infoM$pos[match(pvs$gene, infoM$Composite.Element.REF)]

table(pvs$CN.chr == pvs$gene.chr)
table(pvs$CN.chr == pvs$gene.chr)/nrow(pvs)

dis1 = abs(pvs$gene.pos - pvs$CN.pos) + 1
summary(log10(dis1[which(pvs$CN.chr == pvs$gene.chr)]))

# ------------------------------------------------------------------------
# summary MC relation strength
# ------------------------------------------------------------------------

table(pvs$beta < 0)
table(pvs$beta < 0)/nrow(pvs)
ww1 = which(pvs$CN.chr == pvs$gene.chr)
table(pvs$beta[ww1] < 0)/length(ww1)

ww2 = which(pvs$CN.chr == pvs$gene.chr & dis1 < 1e6)
table(pvs$beta[ww2] < 0)/length(ww2)

setwd("~/research/TCGA/_Sun_MethyE/BRCA/figures2")

pdf("methylation_cn_with_ab_purity_pam50_beta_hist.pdf", width=4, height=3.5)
par(mar=c(5,4,1,1), bty="n")
beta1 = pvs$beta
beta1[which(beta1 > 2)] = 2
beta1[which(beta1 < -2)] = -2
hist(beta1, xlab="beta (methylation vs. SCNA)", main="")
dev.off()

table(pvs$CN.chr)
pvs2 = pvs
pvs2$CN.chr = gsub("chr", "", pvs$CN.chr)
table(pvs2$CN.chr, useNA="ifany")

pvs2$CN.chr[pvs2$CN.chr == "X"] = "23"
pvs2$CN.chr = as.numeric(pvs2$CN.chr)
table(pvs2$CN.chr, useNA="ifany")

pvs2 = pvs2[order(pvs2$CN.chr, pvs2$CN.pos),]
dim(pvs2)
pvs2[1:2,]

ww1 = which(pvs2$CN.chr==1)

pvs2Chr1 = pvs2[ww1,]

n.total = tapply(pvs2Chr1$beta, pvs2Chr1$SNP, length)
n.negtv = tapply(pvs2Chr1$beta, pvs2Chr1$SNP, function(v){length(which(v<0))})

gene1p = infoC$gene[which(infoC$chr=="chr1" & infoC$end < 1.4e8)]
gene1q = infoC$gene[which(infoC$chr=="chr1" & infoC$end >= 1.4e8)]
ww1p   = which(names(n.total) %in% gene1p)
ww1q   = which(names(n.total) %in% gene1q)

xl = "# of CN-methylation association per gene"
yl = "# of negative CN-methy association per gene"

pdf("methylation_cn_with_ab_purity_pam50_chr1.pdf", width=5, height=5)
par(mar=c(5,4,1,1), bty="n")
plot(pvs2$CN.pos[ww1], pvs2$beta[ww1], cex=0.1, xlab="Genoic Location",
    ylab="Beta")

plot(n.total, n.negtv/n.total, cex=0.2, type="n", xlab=xl, ylab=yl)
points(n.total[ww1p], n.negtv[ww1p]/n.total[ww1p], cex=0.3, col="red")
points(n.total[ww1q], n.negtv[ww1q]/n.total[ww1q], cex=0.3, col="blue")
legend("bottomright", legend=c("1p", "1q"), col=c("red", "blue"), pch=1, bty="n")
dev.off()

# ------------------------------------------------------------------------
# check the relation between average copy number changs vs.
# the frequency to see CN-methylation association
# ------------------------------------------------------------------------

n.total = tapply(pvs2$beta, pvs2$SNP, length)
n.negtv = tapply(pvs2$beta, pvs2$SNP, function(v){length(which(v<0))})

pdf("methylation_cn_with_ab_purity_pam50_sign.pdf", width=5, height=5)
par(mar=c(5,4,1,1), bty="n")
plot(n.total, n.negtv/n.total, cex=0.3, xlab=xl, ylab=yl)
abline(h=seq(0,1,by=0.1), col="grey")
abline(v=seq(0,400,by=50), col="grey")
dev.off()

rat = n.negtv/n.total
x2check = which(n.total > 50 & n.total < 150 & rat > 0.5 & rat < 0.7)
g2check = names(n.total)[x2check]

inf2check = infoC[match(g2check, infoC$gene),]
dim(inf2check)
table(inf2check$chr)

# ------------------------------------------------------------------------
# read in CN data
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

datC = read.table("cn_values.txt", sep="\t", header=TRUE, as.is=TRUE)
dim(datC)
datC[1:2,1:5]

ave = apply(datC[,-1], 1, mean)

ave1 = ave[match(names(n.total), datC$id)]

setwd("~/research/TCGA/_Sun_MethyE/BRCA/figures2")

pdf("methylation_cn_with_ab_purity_pam50_ave_CN_vs_sign.pdf", width=8, height=4)

par(mfrow=c(1,2), bty="n", mar=c(5,4,1,1))
plot(ave1, n.total, cex=0.3)
abline(h=seq(0,1,by=0.1), col="grey")
plot(ave1, n.negtv/n.total, cex=0.3)
abline(h=seq(0,1,by=0.1), col="grey")

dev.off()

# ------------------------------------------------------------------------
# whether a methylation site always have positive or negative association
# ------------------------------------------------------------------------

fun1 <- function(v){paste(sort(unique(v)), collapse=";")}
signs = tapply(sign(pvs$beta), pvs$gene, fun1)
length(unique(pvs$gene))
length(signs)
table(signs)

signs[1:5]
mPositive = names(signs)[which(signs=="1")]
mNegative = names(signs)[which(signs=="-1")]

# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# ------------------------------------------------------------------------
# get methylation sites that have postive or negative CN association
# ------------------------------------------------------------------------


table(mPositive %in% infMd$Name)
table(mNegative %in% infMd$Name)

table(mPositive %in% infMd$Name)/length(mPositive)
table(mNegative %in% infMd$Name)/length(mNegative)

mPositive = intersect(mPositive, infMd$Name)
mNegative = intersect(mNegative, infMd$Name)

# ------------------------------------------------------------------------
# are there any significant overlap between CN-assocaited methylation
# sites vs. hot mehyal probes and local methy probes
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

eQTM = read.table("methylation_eQTM_info.txt", sep="\t", header=TRUE,
as.is=TRUE, quote="")
dim(eQTM)
eQTM[1:2,]

length(intersect(eQTM, mPositive))
length(intersect(eQTM, mNegative))

hots = read.table("methylation_hot.txt", sep="\t", header=TRUE,
as.is=TRUE, quote="")
dim(hots)
hots[1:2,]

length(intersect(hots, mPositive))
length(intersect(hots, mNegative))

# ------------------------------------------------------------------------
# generate mehylation probe information
# ------------------------------------------------------------------------

iPositive = infMd[match(mPositive, infMd$Name),]
iNegative = infMd[match(mNegative, infMd$Name),]

dim(iPositive)
dim(iNegative)

iPN = rbind(iPositive, iNegative)
dim(iPN)

methyCNSign = rep(c("+", "-"), times=c(nrow(iPositive), nrow(iNegative)))
methyCNSign[1:5]
table(methyCNSign)

grps = strsplit(iPN$UCSC_RefGene_Group, split=";")
grps = sapply(grps, function(v){paste(sort(unique(v)), collapse=";")})

grps[which(grps=="")] = "None"

t1 = table(grps)
sort(t1, decreasing=TRUE)[1:20]
length(t1)

lbls = c(names(t1)[t1 > 400], "5'UTR;1stExon")
lbls

grps[which(grps=="1stExon;5'UTR")] = "5'UTR;1stExon"
grps[which(! grps %in% lbls)] = "Others"
sort(table(grps))

iPN$Enhancer[which(is.na(iPN$Enhancer))] = "No"
iPN$DHS[which(is.na(iPN$DHS))] = "No"

iPN$Enhancer[which(iPN$Enhancer=="TRUE")] = "Yes"
iPN$DHS[which(iPN$DHS=="TRUE")] = "Yes"

ww1 = which(iPN$Relation_to_UCSC_CpG_Island =="")
iPN$Relation_to_UCSC_CpG_Island[ww1] = "None"

iPN[1:2,]

# ------------------------------------------------------------------------
# check their relation against group
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

chisq.test(methyCNSign, grps)

t1 = table(methyCNSign, grps)
t1
round(t1/rowSums(t1),3)

t1 = t1[,c(7,8,2,3,4,1,6,5)]
t1
nms = colnames(t1)

cols = c("lightcyan", "tomato")

pdf("../figures2/methylation_CN_group.pdf", width=5.5, height=3)
par(mar=c(0.1,2,1,1))
barplot(round(t1/rowSums(t1),3), ylim=c(-0.56,0.52), yaxt="n",
col=cols, beside=TRUE)
axis(side=2, at = seq(0,0.5,by=0.1))
text(x=seq(1.2, length.out=9, by=3), y=rep(-0.02,9), nms, srt=300, pos=4)
text(x=13.2, y=-0.5, "Genomic location group")
legend("topleft", legend=c("mQCN -", "mQCN +"), bty="n", fil=cols)
dev.off()

# ------------------------------------------------------------------------
# check relation against other features
# ------------------------------------------------------------------------

chisq.test(methyCNSign, iPN$Regulatory_Feature_Group)
chisq.test(methyCNSign, iPN$Relation_to_UCSC_CpG_Island)
chisq.test(methyCNSign, iPN$DHS)
chisq.test(methyCNSign, iPN$Enhancer)

t1 = table(methyCNSign, iPN$Regulatory_Feature_Group)
t1
round(t1/rowSums(t1),3)

pdf("../figures2/methylation_CN_CpG.pdf", width=7.5, height=3)
layout(mat=matrix(1:3, nrow=1), widths=c(8,3,3))
par(las=0, mar=c(5,2,1,1), bty="n", cex=1)

t1 = table(methyCNSign, iPN$Relation_to_UCSC_CpG_Island)
t1
t1 = t1[,c(1,3,6,2,5,4)]
t1

colnames(t1)[2:5] = c("North", "South", "North", "South")
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="CpG Island", ylab="Percentage")

t1 = table(methyCNSign, iPN$Enhancer)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="Enhancer", ylab="Percentage")


t1 = table(methyCNSign, iPN$DHS)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="DHS", ylab="Percentage")

dev.off()


q(save = "no")

