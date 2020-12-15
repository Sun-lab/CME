
# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# ------------------------------------------------------------------------
# files
# ------------------------------------------------------------------------

source("~/research/TCGA/_Sun_MethyE/shared_code/eQTL_Plot.R")

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# --------------------------------------------------------------
# read in location annotation
# --------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoE)
infoE[1:2,]

infoM = read.table("methylation_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoM)
infoM[1:2,]

eChr = gsub("chr", "", infoE$chr)
mChr = infoM$Chromosome

table(eChr, useNA="ifany")
table(mChr, useNA="ifany")

cols     = c("green", "blue", "red", "black")
eChr     = gsub("chr", "", infoE$chr)
ePos     = 0.5*(infoE$start + infoE$end)
mChr     = infoM$Chromosome
mPos     = infoM$Genomic_Coordinate
chroms   = 1:22

# --------------------------------------------------------------
# iteratively summarizing results
# --------------------------------------------------------------

inputi = "expression_vs_methylation_rm_cn_PCs7_qnorm.txt"
inputi

outputi = gsub(".txt", ".png", inputi, fixed=TRUE)
outputi

# ------------------------------------------------------------
# read in eQTL results
# ------------------------------------------------------------

pvs = read.table(inputi, header = TRUE, sep = "\t", as.is = TRUE)
dim(pvs)
pvs[1:2,]

table(pvs$SNP  %in% infMd$Name)
pvs = pvs[which(pvs$SNP  %in% infMd$Name),]

summary(pvs$p.value)

peChr = eChr[match(pvs$gene, infoE$gene)]
pmChr = mChr[match(pvs$SNP,  infoM$Composite.Element.REF)]

pcuts = c(0, 10^(seq(-30, -10, by=10)))
np    = length(pcuts) - 1
tbl1  = tbl2 = rep(NA, np)

for (j in 1:np) {
  ww1 = which(pvs$p.value < pcuts[j + 1] & pvs$p.value >= pcuts[j])
  tbl1[j] = length(ww1)
  tbl2[j] = length(which(peChr[ww1] == pmChr[ww1])) / length(ww1)
}
pcuts
cbind(tbl1, tbl2)

# ------------------------------------------------------------
# plot it
# ------------------------------------------------------------

geneID   = match(pvs$gene, infoE$gene)
markerID = match(pvs$SNP,  infoM$Composite.Element.REF)
scores   = pvs$p.value
cols     = c("green", "blue", "red", "black")

# ------------------------------------------------------------
# plot it with more liberal cutoff
# ------------------------------------------------------------

scuts    = c(1e-10, 1e-15, 1e-20, 1e-30)
outputi1 = sub(".png", "_rm_bad_m_qnorm.png", outputi, fixed = TRUE)

png(sprintf("../figures2/extra_liberal_1_%s", outputi1), width = 7.5, 
  height =9, res = 400, units = "in")

eqtl.plot(
  geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab = "Methylation Probe Location",
  ylab = "Transcript Location", plot.hotspots = TRUE,
  hotspots.cut = 10, score.type = "p-value"
)

dev.off()

scuts    = c(1e-20, 1e-25, 1e-30, 1e-40)
outputi1 = sub(".png", "_rm_bad_m_qnorm.png", outputi, fixed = TRUE)

png(sprintf("../figures2/extra_liberal_2_%s", outputi1), width = 7.5, 
    height =9, res = 400, units = "in")

eqtl.plot(
  geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab = "Methylation Probe Location",
  ylab = "Transcript Location", plot.hotspots = TRUE,
  hotspots.cut = 10, score.type = "p-value"
)

dev.off()

q(save = "no")

