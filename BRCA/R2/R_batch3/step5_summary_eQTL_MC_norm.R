
library(data.table)
library(GenomicRanges)

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------------------------
# Read in eQTL results.
# ------------------------------------------------------------------------------

pvs0 = fread("methylation_vs_cn_with_pam50.txt", header=TRUE)
dim(pvs0)
pvs0[1:2,]
summary(pvs0[["p-value"]])

key0 = paste(pvs0$gene, pvs0$SNP, sep=":")

pvs = fread("methylation_vs_cn_with_pam50_qnorm.txt", header=TRUE)
dim(pvs)
pvs[1:2,]
summary(pvs[["p-value"]])

keys = paste(pvs$gene, pvs$SNP, sep=":")

length(intersect(key0, keys))
keyAll = union(key0, keys)
length(keyAll)

p0 = p1 = rep(5, length(keyAll))
p0[match(key0, keyAll)] = -log10(pvs0[["p-value"]])
p1[match(keys, keyAll)] = -log10(pvs[["p-value"]])

png("../figures2/methylation_vs_cn_with_pam50_qnorm_vs_not.png",
      width=5, height=5, res=400, units="in")
par(bty="n", mar=c(5,4,1,1))
plot(p0, p1, xlab="-log10(p-value) Untransformed data",
ylab="-log10(p-value) normalized data", pch=20,
col=rgb(0.8,0.2,0.2,0.5))
abline(0, 1, col="grey", lwd=2)
dev.off()

table(pvs[["p-value"]] < 1e-50)
table(pvs[["p-value"]] < 1e-60)

# ------------------------------------------------------------------------
# files
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

inputs = "methylation_vs_cn_qnorm.txt"
inputs = c(inputs, "methylation_vs_cn_with_pam50_qnorm.txt")
inputs = c(inputs, "methylation_vs_cn_with_ab_purity_pam50_qnorm.txt")
inputs = c(inputs, "methylation_vs_cn_with_pam50_PCs7_qnorm.txt")

inputs

outputs = gsub(".txt", ".png", inputs, fixed=TRUE)

outputs

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

infoE = read.table("methylation_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoE)
infoE[1:2,]

table(infoE$Chromosome, useNA="ifany")

infoM = read.table("cn_info.txt", sep="\t", header=TRUE,
as.is=TRUE)
dim(infoM)
infoM[1:2,]
table(infoM$chr, useNA="ifany")

# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/shared_data/450k")

infMd = read.table("HumanMethylation450_15017482_v1-2_filtered_updated.txt",
sep="\t", header=TRUE, as.is=TRUE, quote="")

dim(infMd)
infMd[1:2,]

# --------------------------------------------------------------
# iteratively summarizing results
# --------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

source("~/research/TCGA/_Sun_MethyE/shared_code/eQTL_Plot.R")

pvsL = list()

for(i in 1:4){
  
  inputi  = inputs[i]
  outputi = outputs[i]
  
  # ------------------------------------------------------------------------
  # read in eQTL results
  # ------------------------------------------------------------------------
  
  pvs = fread(inputi, header=TRUE)
  dim(pvs)
  pvs[1:2,]
  
  summary(pvs[["p-value"]])
  summary(pvs$FDR)
  
  table(pvs[["p-value"]] < 1e-30)
  table(pvs[["p-value"]] < 1e-40)
  table(pvs[["p-value"]] < 1e-50)
  table(pvs[["p-value"]] < 1e-60)
  table(pvs[["p-value"]] < 1e-80)

  table(pvs$gene  %in% infMd$Name)
  pvs = pvs[which(pvs$gene  %in% infMd$Name),]
  pvsL[[i]] = pvs
  
  # ------------------------------------------------------------
  # plot it after removing those possible problematic probes
  # ------------------------------------------------------------
  
  ptplot   = pvs
  
  geneID   = match(ptplot$gene, infoE$Composite.Element.REF)
  markerID = match(ptplot$SNP,  infoM$gene)
  scores   = ptplot[["p-value"]]
  scuts    = c(1e-20, 1e-30, 1e-40, 1e-60)
  cols     = c("green", "blue", "red", "black")
  eChr     = infoE$Chromosome
  ePos     = infoE$Genomic_Coordinate
  mChr     = gsub("chr", "", infoM$chr)
  mPos     = 0.5*(infoM$start + infoM$end)
  chroms   = 1:22
  
  outputi1 = sub(".png", "_rm_bad_m.png", outputi, fixed=TRUE)
  png(sprintf("../figures2/%s", outputi1), width=7.5, height=9,
  res=200, units="in")
  
  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab="Methylation Probe Location",
  ylab="Transcript Location", plot.hotspots=TRUE,
  hotspots.cut=10, score.type="p-value")
  
  dev.off()
}

# --------------------------------------------------------------
# check the intersection
# --------------------------------------------------------------

scuts = c(1e-20, 1e-30, 1e-40, 1e-50, 1e-60)

eMC2 = pvsL[[2]]
eMC3 = pvsL[[3]]

dim(eMC2)
dim(eMC3)
eMC2[1:2,]
eMC3[1:2,]

key2 = paste(eMC2$SNP, eMC2$gene, sep=":")
key3 = paste(eMC3$SNP, eMC3$gene, sep=":")

ns = matrix(nrow=length(scuts), ncol=5)

for(i in 1:length(scuts)){
  s1 = scuts[i]
  
  key2s = key2[which(eMC2[["p-value"]] < s1)]
  key3s = key3[which(eMC3[["p-value"]] < s1)]
  key4s = key3[which(eMC3[["p-value"]] < s1*100)]
  
  ns[i,1] = length(key2s)
  ns[i,2] = length(key3s)
  ns[i,3] = length(key4s)
  
  ns[i,4] = length(intersect(key2s, key3s))
  ns[i,5] = length(intersect(key2s, key4s))
}

ns
ns/ns[,1]

# --------------------------------------------------------------
# check the hot spot at chromosome 16
# --------------------------------------------------------------

pvs = pvsL[[2]]
dim(pvs)
pvs[1:5,]

cpgID     = match(pvs$gene, infoE$Composite.Element.REF)
geneID    = match(pvs$SNP,  infoM$gene)
pvs$mChr  = infoE$Chromosome[cpgID]
pvs$mPos  = infoE$Genomic_Coordinate[cpgID]
pvs$mGene = infoE$Gene_Symbol[cpgID]
pvs$eChr  = gsub("chr", "", infoM$chr)[geneID]
pvs$ePos  = 0.5*(infoM$start + infoM$end)[geneID]

dim(pvs)
pvs[1:5,]

w2check = which(pvs$eChr=="17" & pvs$mChr != "17")
length(w2check)

w2check = which(pvs$eChr=="16" & pvs$mChr != "16")
length(w2check)

pvs2chk = pvs[w2check,]
dim(pvs2chk)
pvs2chk[1:50,-(3:5)]
summary(pvs2chk$beta)

sort(table(pvs2chk$mGene), decreasing=TRUE)[1:50]
sort(table(pvs2chk$SNP), decreasing=TRUE)[1:50]

tt1 = sort(table(pvs2chk$SNP), decreasing=TRUE)
gene2check = names(tt1)[which(tt1 > 90)]
length(gene2check)

geneSyms   = infoM$geneSymbol[match(gene2check, infoM$gene)]
cat(geneSyms, sep="\n")

geneInfo = infoM[match(gene2check,infoM$gene),-c(1, 8)]
dim(geneInfo)
geneInfo[1:2,]

geneInfo$freq = as.numeric(tt1[match(gene2check, names(tt1))])
dim(geneInfo)
geneInfo[1:2,]

write.table(geneInfo, sep=" & ", row.names = FALSE, eol = "\\\\\n", 
            quote=FALSE)

# ------------------------------------------------------------------------------
# Read in results of SCNA vs. gene expression associations
# ------------------------------------------------------------------------------

pvsEC = fread("expression_vs_cn_with_pam50_qnorm.txt", header=TRUE)
dim(pvsEC)
pvsEC[1:2,]

table(pvsEC$SNP == pvsEC$gene)
pvsEC = pvsEC[which(pvsEC$SNP == pvsEC$gene),]
dim(pvsEC)
pvsEC[1:2,]

summary(pvsEC[["p-value"]])

eID       = match(pvsEC$SNP,  infoM$gene)
pvsEC$Chr  = gsub("chr", "", infoM$chr)[eID]
pvsEC$Pos  = 0.5*(infoM$start + infoM$end)[eID]
dim(pvsEC)
pvsEC[1:2,]

wEC = which(pvsEC$Chr=="16")
length(wEC)

pvsEC2chk = pvsEC[wEC,]
dim(pvsEC2chk)

table(gene2check %in% pvsEC2chk$gene)
pvsEC2chk[which(pvsEC2chk$gene %in% gene2check),]

tb2 = pvsEC2chk[which(pvsEC2chk$gene %in% gene2check),]
mat1 = match(tb2$SNP, infoM$gene)
tb3  = cbind(infoM[mat1,c(2,4:7)], signif(tb2[[5]],2))
dim(tb3)
tb3[1:5,]

names(tb3)[6] = "p-value"

write.table(tb3, sep=" & ", row.names = FALSE, eol = "\\\\\n", 
            quote=FALSE)

dim(tb3)

# ------------------------------------------------------------------------------
# Read in SCNA data and methylation data.
# ------------------------------------------------------------------------------

datC = read.table(file = "cn_values.txt", sep = "\t",
                  header = TRUE, as.is = TRUE)
dim(datC)
datC[1:2, 1:5]

datM = fread(file = "methylation_mvalue_qnorm.txt")
dim(datM)
datM[1:2, 1:5]

table(names(datC) == names(datM))

# datM is too big to do analysis below; restrict to only the relevant probes.
probes_to_keep = unique(pvs2chk$gene)

dim(pvs2chk)
length(probes_to_keep)

dim(datM)
datM = datM[which(datM$id %in% probes_to_keep), ]
dim(datM)

pvs2chk = pvs2chk[order(pvs2chk[["p-value"]]),]
dim(pvs2chk)
pvs2chk[1:5,]

cn1 = as.numeric(datC[which(datC$id == "RANBP10|57610"),-1])
me1 = as.numeric(datM[which(datM$id == "cg22385764"),-1])

# ------------------------------------------------------------------------------
# check those CpG's
# ------------------------------------------------------------------------------

dim(pvs2chk)
pvs2chk[1:2,]
summary(pvs2chk$beta)
summary(pvs2chk[["p-value"]])

cpgs = unique(pvs2chk$gene)
length(cpgs)

genes = unique(pvs2chk$SNP)
length(genes)

t.cpgs  = table(pvs2chk$gene)
t.genes = table(pvs2chk$SNP)
summary(as.numeric(t.cpgs))
summary(as.numeric(t.genes))

dim(infMd)
infMd[1:2,]

table(cpgs %in% infMd$Name)
wCpg = which(infMd$Name %in% cpgs)
length(wCpg)

table(infMd$DHS, useNA="ifany")
table(infMd$DHS[wCpg], useNA="ifany")

cpg16 = infMd$Name %in% cpgs
c1 = chisq.test(!is.na(infMd$DHS), cpg16)
c1
round(c1$expected)
c1$observed

table(infMd$Relation_to_UCSC_CpG_Island, useNA="ifany")
table(infMd$Relation_to_UCSC_CpG_Island[wCpg], useNA="ifany")

c2 = chisq.test(infMd$Relation_to_UCSC_CpG_Island, cpg16)
c2
round(c2$expected)
c2$observed

# ------------------------------------------------------------------------------
# read in CTCF binding sites
# ------------------------------------------------------------------------------

info.cpg16 = infMd[cpg16, c("Name", "CHR", "MAPINFO", "MAPINFO")]
dim(info.cpg16)
info.cpg16[1:2,]
names(info.cpg16) = c("name", "chr", "start", "end")
info.cpg16$chr = paste0("chr", info.cpg16$chr)

dim(info.cpg16)
info.cpg16[1:2,]

gr1 = makeGRangesFromDataFrame(info.cpg16, ignore.strand=TRUE)
gr1

fun1 <- function(x) sum(width(reduce(x, ignore.strand=T)))
fun1(gr1)

setwd("~/research/data/CTCF/")

ff1  = "CTCFBSDB_all_exp_sites_Sept12_2012_hg19_loci.bed"
ctcf = fread(ff1)
dim(ctcf)
ctcf[1:5,]

table(ctcf$chr)
table(info.cpg16$chr)

names(ctcf) = c("chr", "start", "end")

lens = ctcf$end - ctcf$start + 1
pdf("CTCFBS_len_hist.pdf", width=6, height=4)
par(mar=c(5,4,1,1), bty="n")
hist(log10(lens), xlab="log10(CTCF BS length)", main="", breaks=100)
abline(v=log10(350))
dev.off()

sum(lens)
table(lens < 500)
table(lens < 400)
table(lens < 350)

sum(lens)
sum(lens[which(lens < 350)])

gr2 = makeGRangesFromDataFrame(ctcf, ignore.strand=TRUE, seqnames.field="chr",
                               start.field="start", end.field="end")

gr3 = makeGRangesFromDataFrame(ctcf[which(lens < 350),], ignore.strand=TRUE, 
                               seqnames.field="chr",
                               start.field="start", end.field="end")

width2 = fun1(gr2)
width3 = fun1(gr3)
width2
width3

prp2 = width2/(3234.83*10^6)
prp3 = width3/(3234.83*10^6)
prp2
prp3

mtch2 = findOverlaps(gr1, gr2, select="first")
table(!is.na(mtch2))

mtch3 = findOverlaps(gr1, gr3, select="first")
table(!is.na(mtch3))

pbinom(134, 123-1, prp2, lower.tail=FALSE, log.p=TRUE)
pbinom(134, 120-1, prp3, lower.tail=FALSE, log.p=TRUE)

q(save="no")

