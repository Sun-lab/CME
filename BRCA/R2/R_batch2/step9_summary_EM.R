
# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

info = read.table("methylation_eQTM_info.txt", sep="\t", header=TRUE,
                    as.is=TRUE, quote="")
dim(info)
info[1:2,]

pval = read.table("methylation_eQTM_pval.txt", sep="\t", header=TRUE,
                    as.is=TRUE, quote="")

dim(pval)
pval[1:2,]

table(pval$SNP %in% info$Name)
pval = pval[which(pval$SNP %in% info$Name),]
dim(pval)

table(info$Genome_Build)

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
                   as.is=TRUE)
dim(infoE)
infoE[1:2,]

table(infoE$chr, useNA="ifany")
infoE$chr = gsub("chr", "", infoE$chr)
table(infoE$chr, useNA="ifany")

infoM = read.table("methylation_info.txt", sep="\t", header=TRUE,
                   as.is=TRUE)
dim(infoM)
infoM[1:2,]
table(infoM$Chromosome, useNA="ifany")

# ------------------------------------------------------------------------
# read in gene location information from emsembl
# ------------------------------------------------------------------------

ensembl = read.table("~/research/data/human/hg19_info/ensGene_info.bed",
  sep="\t", as.is=TRUE)

dim(ensembl)
ensembl[1:5,]

table(is.na(infoE$ensembl))
enID = intersect(infoE$ensembl, ensembl$V4)
length(enID)

infoE1 = infoE[match(enID, infoE$ensembl),]
infoE2 = ensembl[match(enID, ensembl$V4),]

dim(infoE1)
dim(infoE2)
head(infoE1)
head(infoE2)

table(infoE2$V2 < infoE2$V3)
table(paste0("chr", infoE1$chr) == infoE2$V1)

table(infoE1$strand == infoE2$V6)
table(infoE1$start - infoE2$V2 > 0)
table(infoE1$end - infoE2$V3 > 0)
summary(infoE1$start - infoE2$V2)
summary(infoE1$end - infoE2$V3)

# the location information of array annotation and ensembl are
# very similar but not the same. We will still use the array annotation
# because the array probe may only detect some transcripts of a gene

# ------------------------------------------------------------------------
# check the distance between methy probe and associated genes
# ------------------------------------------------------------------------

table(pval$SNP  %in% infoM$Composite.Element.REF)
table(pval$SNP  %in% info$Name)
table(pval$gene %in% infoE$gene)

iM0 = infoM[match(pval$SNP,  infoM$Composite.Element.REF),]
iM  = info[match(pval$SNP,   info$Name),]
iE  = infoE[match(pval$gene, infoE$gene),]

iM0[1:5,]
iM[1:5,]
iE[1:5,]

table(iM0$Genomic_Coordinate == iM$MAPINFO)
table(iM$CHR == iM0$Chromosome)

table(iM$CHR == iE$chr)

dis = rep(NA, nrow(iM))
ww1 = which(iM$CHR == iE$chr)

summary(iE$end - iE$start)
table(iE$strand, useNA="ifany")

TSS        = iE$start
wnegS      = which(iE$strand == "-")
TSS[wnegS] = iE$end[wnegS]

dis[ww1]   = iM$MAPINFO[ww1] - TSS[ww1]
dis[wnegS] = -dis[wnegS]

# ------------------------------------------------------------------------
# standardize distance by gene length
# ------------------------------------------------------------------------

disNew  = dis
geneLen = iE$end - iE$start

wpp = which(dis>0)
disNew[wpp] = dis[wpp]/geneLen[wpp]

wnn = which(dis<0)
disNew[wnn] = dis[wnn]/1000

quantile(disNew[wpp], probs = seq(0, 1, 0.05))
quantile(disNew[wnn], probs = seq(0, 1, 0.05))

table(disNew < -1)
table(disNew > 1)
w2kp   = which(disNew <= 1 & disNew >= -1)
length(w2kp)

disClass = rep(NA, length(disNew))
disClass[which(disNew < -1)] = "TSS1k+"
disClass[which(disNew <= 0 & disNew >= -1)] = "TSS1k"
disClass[which(disNew > 0 & disNew <= 1)] = "Body"
disClass[which(disNew > 1)] = "3End"

signs = rep("-", nrow(pval))
signs[which(pval$beta > 0)] = "+"

t1 = table(signs, disClass)
t1
t1 = t1[,4:1]
round(t1/rowSums(t1),2)

cols = c("palegreen", "seashell")

pdf("../figures2/methylation_eQTM_to_TSS_barplot.pdf", width=4, height=3)
par(mar=c(2.5,4,1,1))
barplot(round(t1/rowSums(t1),2), beside=T, col=cols, ylab="Proportion")
legend("topleft", legend=c("eQTM -", "eQTM +"), fill=cols, bty="n")
dev.off()

wPositive = which(pval$beta[w2kp] > 0)
wNegative = which(pval$beta[w2kp] < 0)
propPositive = length(wPositive)/length(w2kp)
propPositive

d0 = density(disNew[w2kp])
d1 = density(disNew[w2kp][wPositive])
d2 = density(disNew[w2kp][wNegative])

pdf("../figures2/methylation_eQTM_to_TSS.pdf", width=6, height=4)
par(mar=c(5,4,0,1))
plot(d2$x, d2$y, xlab="(scaled) Distance to TSS", main="", bty="n",
col="darkblue", lwd=2, lty=1, type="l", ylab="Density")
lines(d1$x, d1$y, col="darkred", lwd=2, lty=2)
abline(v=0)

legend("topright", legend=c("eQTM +", "eQTM -"), bty="n",
lty=c(2,1), lwd=c(2,2), col=c("darkred", "darkblue"))

dev.off()

# ------------------------------------------------------------------------
# check the annotation of those with positve and negative association
# ------------------------------------------------------------------------

chisq.test(pval$beta > 0, iM$CHR == iE$chr)
t1 = table(pval$beta > 0, iM$CHR == iE$chr)
t1
t1/rowSums(t1)

methyPo = unique(pval$SNP[pval$beta > 0])
methyNe = unique(pval$SNP[pval$beta < 0])

length(methyPo)
length(methyNe)

methBo = intersect(methyPo, methyNe)
length(methBo)

methyPo = setdiff(methyPo, methBo)
methyNe = setdiff(methyNe, methBo)

length(methyPo)
length(methyNe)

info1 = info[-which(info$Name %in% methBo),]
dim(info1)
info1[1:2,]

sign = rep("+", nrow(info1))
sign[which(info1$Name %in% methyNe)] = "-"
table(info1$Name[which(sign=="+")] %in% methyPo)

# ------------------------------------------------------------------------
# perform test
# ------------------------------------------------------------------------

grps = strsplit(info1$UCSC_RefGene_Group, split=";")
grps = sapply(grps, function(v){paste(sort(unique(v)), collapse=";")})

grps[which(grps=="")] = "None"

t1 = table(grps)
sort(t1, decreasing=TRUE)[1:20]
length(t1)

lbls = c(names(t1)[t1 > 200], "5'UTR;1stExon", "3'UTR")
lbls

grps[which(grps=="1stExon;5'UTR")] = "5'UTR;1stExon"
grps[which(! grps %in% lbls)] = "Others"
sort(table(grps))

# ------------------------------------------------------------------------
# check their relation against group
# ------------------------------------------------------------------------

chisq.test(sign, grps)

t1 = table(sign, grps)
t1

t1 = t1[,c(7,8,2,3,4,1,6,5)]
t1
round(t1/rowSums(t1),3)

nms = colnames(t1)

cols = c("palegreen", "seashell")

pdf("../figures2/methylation_eQTM_sign_group.pdf", width=5.5, height=3)
par(mar=c(0.1,2,1,1))
barplot(round(t1/rowSums(t1),3), ylim=c(-0.55,0.6), yaxt="n",
col=cols, beside=TRUE)
axis(side=2, at = seq(0,0.6,by=0.1))
text(x=seq(1.2, length.out=9, by=3), y=rep(-0.02,9), nms, srt=300, pos=4)
text(x=13.2, y=-0.5, "Genomic location group")
legend("topleft", legend=c("eQTM -", "eQTM +"), bty="n", fil=cols)
dev.off()

# ------------------------------------------------------------------------
# check their relation against other features
# ------------------------------------------------------------------------

chisq.test(sign, info1$Regulatory_Feature_Group)
chisq.test(sign, info1$Relation_to_UCSC_CpG_Island)

chisq.test(sign, info1$Enhancer)
chisq.test(sign, info1$DHS)

t1 = table(sign, info1$Regulatory_Feature_Group)
t1
round(t1/rowSums(t1),3)

pdf("../figures2/methylation_eQTM_sign_CpG.pdf", width=7.5, height=3)
layout(mat=matrix(1:3, nrow=1), widths=c(8,3,3))
par(las=0, mar=c(5,2,1,1), bty="n", cex=1)

t1 = table(sign, info1$Relation_to_UCSC_CpG_Island)
t1
t1 = t1[,c(1,3,6,2,5,4)]
colnames(t1)[2:5] = c("North", "South", "North", "South")
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="CpG Island", ylab="Percentage")

t1 = table(sign, info1$Enhancer)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="Enhancer", ylab="Percentage")


t1 = table(sign, info1$DHS)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="DHS", ylab="Percentage")

dev.off()


q(save = "no")

