
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoE)
infoE[1:2,]
infoE$pos = 0.5*(infoE$start + infoE$end)

table(infoE$chr, useNA="ifany")

infoC = read.table("cn_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoC)
infoC[1:2,]
infoC$pos = 0.5*(infoC$start + infoC$end)

table(infoC$chr, useNA="ifany")

# ------------------------------------------------------------------------
# read in association results
# ------------------------------------------------------------------------

pvs0 = read.table("expression_vs_cn_qnorm.txt",
                  header=TRUE, sep="\t", as.is=TRUE)
dim(pvs0)
pvs0[1:2,]

summary(pvs0$p.value)
summary(pvs0$FDR)

table(pvs0$p.value < 1e-30)

pvs0$CN.chr = infoC$chr[match(pvs0$SNP, infoC$gene)]
pvs0$CN.pos = infoC$pos[match(pvs0$SNP, infoC$gene)]

pvs0$gene.chr = infoE$chr[match(pvs0$gene, infoE$gene)]
pvs0$gene.pos = infoE$pos[match(pvs0$gene, infoE$gene)]

table(pvs0$CN.chr == pvs0$gene.chr)
table(pvs0$CN.chr == pvs0$gene.chr)/nrow(pvs0)

# ------------------------------------------------------------------------
# read in eQTL results, given purity and subtype
# ------------------------------------------------------------------------

pvs = read.table("expression_vs_cn_with_ab_purity_pam50_qnorm.txt",
header=TRUE, sep="\t", as.is=TRUE)
dim(pvs)
pvs[1:2,]

summary(pvs$p.value)
summary(pvs$FDR)

table(pvs$p.value < 1e-30)

pvs$CN.chr = infoC$chr[match(pvs$SNP, infoC$gene)]
pvs$CN.pos = infoC$pos[match(pvs$SNP, infoC$gene)]

pvs$gene.chr = infoE$chr[match(pvs$gene, infoE$gene)]
pvs$gene.pos = infoE$pos[match(pvs$gene, infoE$gene)]

table(pvs$CN.chr == pvs$gene.chr)
table(pvs$CN.chr == pvs$gene.chr)/nrow(pvs)

table(pvs$CN.chr == pvs$gene.chr, pvs$p.value < 1e-30)
table(pvs$CN.chr == pvs$gene.chr, pvs$p.value < 1e-50)

# ------------------------------------------------------------------------
# comparison
# ------------------------------------------------------------------------

summary(pvs0$beta)
table(pvs0$beta > 0)

summary(pvs$beta)
table(pvs$beta > 0)

pvs2ck = pvs[which(pvs$beta < 0),]
pvs2ck[order(pvs2ck$CN.chr, pvs2ck$CN.pos),]

q(save = "no")

