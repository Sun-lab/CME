
# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

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

# ------------------------------------------------------------
# compare old and new methylation information
# ------------------------------------------------------------

table(infoM$Composite.Element.REF %in% infMd$Name)

methyProbe = intersect(infoM$Composite.Element.REF, infMd$Name)
methInfo0  = infoM[match(methyProbe, infoM$Composite.Element.REF),]
methInfo1  = infMd[match(methyProbe, infMd$Name),]

dim(methInfo0)
dim(methInfo1)

methInfo0[1:2,]
methInfo1[1:2,]

table(methInfo0$Composite.Element.REF == methInfo1$Name)

table(methInfo0$Chromosome == methInfo1$CHR)
table(methInfo0$Genomic_Coordinate == methInfo1$MAPINFO)
table(methInfo0$Gene_Symbol == methInfo1$UCSC_RefGene_Name)

ww1 = which(methInfo0$Gene_Symbol != methInfo1$UCSC_RefGene_Name)
methInfo0[ww1[1:5],]
methInfo1[ww1[1:5],]

# ------------------------------------------------------------
# update methylation information
# ------------------------------------------------------------

grps = methInfo1$UCSC_RefGene_Group
grps = strsplit(methInfo1$UCSC_RefGene_Group, split=";")
grps = sapply(grps, function(v){paste(sort(unique(v)), collapse=";")})

grps[which(grps=="")] = "None"

t1 = table(grps)
sort(t1, decreasing=TRUE)[1:20]
length(t1)

lbls = c(names(t1)[t1 > 10000], "5'UTR;1stExon")
lbls

grps[which(grps=="1stExon;5'UTR")] = "5'UTR;1stExon"
grps[which(! grps %in% lbls)] = "Others"
sort(table(grps))

methInfo1$Enhancer[which(is.na(methInfo1$Enhancer))] = "No"
methInfo1$DHS[which(is.na(methInfo1$DHS))] = "No"

methInfo1$Enhancer[which(methInfo1$Enhancer=="TRUE")] = "Yes"
methInfo1$DHS[which(methInfo1$DHS=="TRUE")] = "Yes"

ww1 = which(methInfo1$Relation_to_UCSC_CpG_Island =="")
methInfo1$Relation_to_UCSC_CpG_Island[ww1] = "None"

methInfo1[1:2,]

# ------------------------------------------------------------------------
# read in results summary
# ------------------------------------------------------------------------

gene1e60 = read.table("gene1e60.txt", header=TRUE, sep="\t", as.is=TRUE)
meth1e60 = read.table("meth1e60.txt", header=TRUE, sep="\t", as.is=TRUE)

dim(gene1e60)
gene1e60[c(1:5,96:100,360:362),]

dim(meth1e60)
meth1e60[c(1:5,nrow(meth1e60)),]

summary(gene1e60$freqIn1e60)
summary(meth1e60$freqIn1e60)

table(gene1e60$freqIn1e60 == sort(gene1e60$freqIn1e60, decreasing=TRUE))
table(meth1e60$freqIn1e60 == sort(meth1e60$freqIn1e60, decreasing=TRUE))

# ------------------------------------------------------------------------
# select genes that appear at least 100 times in the 1e60 list
# ------------------------------------------------------------------------

summary(gene1e60$freqIn1e60)
gene100 = gene1e60[gene1e60$freqIn1e60 >= 100,]
dim(gene100)
gene100[c(1:5,nrow(gene100)),]

table(gene100$gene %in% infoE$gene)

gene100 = cbind(gene100, infoE[match(gene100$gene, infoE$gene),-1])
dim(gene100)
gene100[1:5,]

table(gene100$chr)
table(is.na(gene100$ensembl))

cat(gene100$ensembl, sep="\n")

# ------------------------------------------------------------------------
# read in gene expression data
# ------------------------------------------------------------------------

eDat = read.table("expression_log_TReC_rm_cn.txt", header=TRUE,
                  sep="\t", as.is=TRUE)
dim(eDat)
eDat[1:2,1:5]

xDat = read.table("cov_EM_with_ab_purity_pam50.txt", header=TRUE,
as.is=TRUE, sep="\t")

dim(xDat)
xDat[1:2,1:5]

table(names(eDat) == names(xDat))

ids  = xDat$id
xDat = t(xDat[,-1])
xDat = data.frame(xDat)
names(xDat) = ids
dim(xDat)
xDat[1:2,1:5]

names(xDat)

# ------------------------------------------------------------------------
# check the assocaiton of these genes with purity
# ------------------------------------------------------------------------

ePdat = eDat[match(gene100$gene, eDat$id),]
dim(ePdat)
ePdat[1:2,1:5]

table(gene100$gene == ePdat$id)
ePdat = data.matrix(ePdat[,-1])

beta = pvs = rep(NA, nrow(gene100))

for(i in 1:nrow(gene100)){
  gi = gene100$gene[i]
  yi = ePdat[i,]
  
  li = lm(yi ~., data=xDat)
  si = summary(li)
  wi = which(rownames(si$coef) == "purity")
  beta[i] = si$coef[wi,1]
  pvs[i]  = si$coef[wi,4]
}

summary(beta)
summary(-log10(pvs))

pdf("../figures2/hot_190_genes_vs_purity.pdf", width=4, height=4)
par(mar=c(5,4,1,1), bty="n")
plot(beta, -log10(pvs), ylab="-log10(pvalues)", pch=20, col=rgb(1,0.1,0.1,0.6))
dev.off()

# ------------------------------------------------------------------------
# select top methyaltion sites that appear most frequently
# in the 1e60 list
# ------------------------------------------------------------------------

summary(meth1e60$freqIn1e60)
table(meth1e60$freqIn1e60 >= 100)
table(meth1e60$freqIn1e60 >= 50)
table(meth1e60$freqIn1e60 >= 30)
table(meth1e60$freqIn1e60 >= 20)

meth30 = meth1e60[which(meth1e60$freqIn1e60 >= 30),]
dim(meth30)
meth30[c(1:5,nrow(meth30)),]

table(meth30$methyProbe %in% infoM$Composite.Element.REF)
table(meth30$methyProbe %in% infMd$Name)

meth30 = meth30[which(meth30$methyProbe %in% methInfo1$Name),]
dim(meth30)

mat1 = match(meth30$methyProbe, methInfo1$Name)
meth30 = cbind(meth30, methInfo1[mat1,-1])
dim(meth30)
meth30[1:5,]

table(meth30$CHR, useNA="ifany")

hotMeth = rep("No", nrow(methInfo1))
hotMeth[which(methInfo1$Name %in% meth30$methyProbe)] = "Yes"
table(hotMeth)

write.table(meth30, file = "methylation_hot.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = FALSE,
  col.names = TRUE)

# ------------------------------------------------------------------------
# check their relation against group
# ------------------------------------------------------------------------

chisq.test(hotMeth, grps)

t1 = table(hotMeth, grps)
t1
round(t1/rowSums(t1),3)

t1 = t1[,c(7,8,2,3,4,1,6,5)]
t1
nms = colnames(t1)

cols = c("lightskyblue", "honeydew")

pdf("../figures2/methylation_hot_group.pdf", width=5.5, height=3)
par(mar=c(0.1,2,1,1))
barplot(round(t1/rowSums(t1),3), ylim=c(-0.36,0.35), yaxt="n",
    col=cols, beside=TRUE)
axis(side=2, at = seq(0,0.3,by=0.1))
text(x=seq(1.2, length.out=9, by=3), y=rep(-0.02,9), nms, srt=300, pos=4)
text(x=13.2, y=-0.33, "Genomic location group")
legend("topleft", legend=c("Cold", "Hot"), bty="n", fil=cols)
dev.off()

# ------------------------------------------------------------------------
# check relation against other features
# ------------------------------------------------------------------------

chisq.test(hotMeth, grps)
chisq.test(hotMeth, methInfo1$Enhancer)
chisq.test(hotMeth, methInfo1$Regulatory_Feature_Group)
chisq.test(hotMeth, methInfo1$DHS)

t1 = table(hotMeth, methInfo1$Regulatory_Feature_Group)
t1
round(t1/rowSums(t1),3)

pdf("../figures2/methylation_hot_CpG.pdf", width=7.5, height=3)
layout(mat=matrix(1:3, nrow=1), widths=c(8,3,3))
par(las=0, mar=c(5,2,1,1), bty="n", cex=1)

t1 = table(hotMeth, methInfo1$Relation_to_UCSC_CpG_Island)
t1
t1 = t1[,c(1,3,6,2,5,4)]
colnames(t1)[2:5] = c("North", "South", "North", "South")
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
      xlab="CpG Island", ylab="Percentage")

t1 = table(hotMeth, methInfo1$Enhancer)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
        xlab="Enhancer", ylab="Percentage")


t1 = table(hotMeth, methInfo1$DHS)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="DHS", ylab="Percentage")

dev.off()

# ------------------------------------------------------------------------
# read eQTM results
# ------------------------------------------------------------------------

pvs3 = read.table("expression_vs_methylation_ECM_purity7.txt",
header=TRUE, sep="\t", as.is=TRUE)
dim(pvs3)
pvs3[1:2,]

summary(pvs3$p.value)
table(pvs3$p.value < 1e-30)

iE = infoE[match(pvs3$gene, infoE$gene),]
iM = infoM[match(pvs3$SNP,  infoM$Composite.Element.REF),]

dim(iE)
dim(iM)

iE[1:2,]
iM[1:2,]

summary(pvs3$beta)
summary(pvs3$beta[pvs3$p.value < 1e-30])

table(pvs3$beta < 0)
table(pvs3$beta[pvs3$p.value < 1e-30] < 0)

brks = seq(-2, 2, by=0.05)

pdf("../figures2/expression_vs_methylation_beta.pdf", width=7, heigh=6)
par(mfrow=c(2,1), mar=c(5,4,2,1))

hist(pvs3$beta[abs(pvs3$beta)<2], xlab="beta", col=rgb(0.8, 0.1, 0.1, 0.5),
      main="p-value < 1e-10", breaks=brks, border=rgb(0.8, 0.1, 0.1))

legend("topright", bty="n", legend=c("all eQTM", "same-chromosome eQTM"),
fill=c(rgb(0.8, 0.1, 0.1, 0.5), rgb(0.1, 0.1, 0.8, 0.5)),
border=c(rgb(0.8, 0.1, 0.1), rgb(0.1, 0.1, 0.8)))


hist(pvs3$beta[which(abs(pvs3$beta)<2 & iE$chr == iM$Chromosome)],
      xlab="beta", main="p-value < 1e-30 & same chromosome", add=TRUE,
      breaks=brks, col=rgb(0.1, 0.1, 0.8, 0.5), border=rgb(0.1, 0.1, 0.8))

hist(pvs3$beta[which(pvs3$p.value < 1e-30)], xlab="beta",
      col=rgb(0.8, 0.1, 0.1, 0.5), main="p-value < 1e-30", breaks=brks,
      border=rgb(0.8, 0.1, 0.1))

hist(pvs3$beta[which(pvs3$p.value < 1e-30 & iE$chr == iM$Chromosome)],
      xlab="beta", main="p-value < 1e-30 & same chromosome", add=TRUE,
      breaks=brks, col=rgb(0.1, 0.1, 0.8, 0.5),
      border=rgb(0.1, 0.1, 0.8))

dev.off()

# ------------------------------------------------------------------------
# methylation probes with EM association after EMC purity correction
# ------------------------------------------------------------------------

pv4 = pvs3[which(pvs3$p.value < 1e-30),]
dim(pv4)
pv4[1:2,]

length(unique(pv4$SNP))
length(unique(pv4$gene))

sort(table(pv4$SNP), decreasing=TRUE)[1:20]
sort(table(pv4$gene), decreasing=TRUE)[1:20]

table(unique(pv4$SNP) %in% methInfo1$Name)

methLocal = methInfo1[which(methInfo1$Name %in% unique(pv4$SNP)),]
dim(methLocal)
methLocal[1:5,]

table(methLocal$CHR, useNA="ifany")

localMeth = rep("No", nrow(methInfo1))
localMeth[which(methInfo1$Name %in% unique(pv4$SNP))] = "Yes"
table(localMeth)

chisq.test(localMeth, hotMeth)
table(localMeth, hotMeth)

write.table(methLocal, file = "methylation_eQTM_info.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

write.table(pv4, file = "methylation_eQTM_pval.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------------------
# check their relation against group
# ------------------------------------------------------------------------

chisq.test(localMeth, grps)

t1 = table(localMeth, grps)
t1
round(t1/rowSums(t1),3)

t1 = t1[,c(7,8,2,3,4,1,6,5)]
t1
nms = colnames(t1)

cols = c("lightskyblue", "honeydew")

pdf("../figures2/methylation_eQTM_group.pdf", width=5.5, height=3)
par(mar=c(0.1,2,1,1))
barplot(round(t1/rowSums(t1),3), ylim=c(-0.36,0.35), yaxt="n",
col=cols, beside=TRUE)
axis(side=2, at = seq(0,0.3,by=0.1))
text(x=seq(1.2, length.out=9, by=3), y=rep(-0.02,9), nms, srt=300, pos=4)
text(x=13.2, y=-0.33, "Genomic location group")
legend("topleft", legend=c("None eQTM", "eQTM"), bty="n", fil=cols)
dev.off()

# ------------------------------------------------------------------------
# check relation against other features
# ------------------------------------------------------------------------

chisq.test(localMeth, methInfo1$Relation_to_UCSC_CpG_Island)
chisq.test(localMeth, methInfo1$Enhancer)
chisq.test(localMeth, methInfo1$Regulatory_Feature_Group)
chisq.test(localMeth, methInfo1$DHS)

t1 = table(localMeth, methInfo1$Regulatory_Feature_Group)
t1
round(t1/rowSums(t1),3)


pdf("../figures2/methylation_eQTM_CpG.pdf", width=7.5, height=3)
layout(mat=matrix(1:3, nrow=1), widths=c(8,3,3))
par(las=0, mar=c(5,2,1,1), bty="n", cex=1)

t1 = table(localMeth, methInfo1$Relation_to_UCSC_CpG_Island)
t1
t1 = t1[,c(1,3,6,2,5,4)]
colnames(t1)[2:5] = c("North", "South", "North", "South")
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="CpG Island", ylab="Percentage")

t1 = table(localMeth, methInfo1$Enhancer)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="Enhancer", ylab="Percentage")

t1 = table(localMeth, methInfo1$DHS)
t1
round(t1/rowSums(t1),3)

barplot(round(t1/rowSums(t1),3), beside=TRUE, col=cols,
xlab="DHS", ylab="Percentage")

dev.off()

q(save = "no")

