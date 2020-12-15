
# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------------------
# read in results summary
# ------------------------------------------------------------------------

gene1e60 = read.table("gene1e60_qnorm.txt", header=TRUE, sep="\t", as.is=TRUE)
meth1e60 = read.table("meth1e60_qnorm.txt", header=TRUE, sep="\t", as.is=TRUE)

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

meth30 = meth1e60[which(meth1e60$freqIn1e60 >= 30),]
dim(meth30)
meth30[c(1:5,nrow(meth30)),]

# ------------------------------------------------------------------------
# read in gene expression data
# ------------------------------------------------------------------------

eDat = read.table("expression_log_TReC.txt", header=TRUE,
sep="\t", as.is=TRUE)
dim(eDat)
eDat[1:2,1:5]

table(gene100$gene %in% eDat$id)

gd = data.matrix(eDat[match(gene100$gene, eDat$id),-1])
dim(gd)
gd[1:2,1:5]

pdf("../figures2/hot_genes_expression.pdf", width=8, height=4)
par(mfrow=c(1,2))
for(i in 1:nrow(gd)){
  hist(gd[i,], xlab="log(TReC)", main="")
  hist(exp(gd[i,]), xlab="TReC", main="")
}
dev.off()

# ------------------------------------------------------------------------
# read in DNA methylation data
# ------------------------------------------------------------------------

mDat = read.table("methylation_mvalue.txt", header=TRUE,
sep="\t", as.is=TRUE)
dim(mDat)
mDat[1:2,1:5]

table(meth30$methyProbe %in% mDat$id)

md = data.matrix(mDat[match(meth30$methyProbe, mDat$id),-1])
dim(md)
md[1:2,1:5]

beta = exp(md)/(1 + exp(md))

pdf("../figures2/hot_methy.pdf", width=8, height=4)
par(mfrow=c(1,2))
for(i in 1:nrow(gd)){
  hist(md[i,], xlab="M-value", main="")
  hist(beta[i,], xlab="beta-value", main="")
}
dev.off()



q(save = "no")

