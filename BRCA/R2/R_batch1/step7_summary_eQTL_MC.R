
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

inputs = "methylation_vs_cn.txt"
inputs = c(inputs, "methylation_vs_cn_with_pam50.txt")
inputs = c(inputs, "methylation_vs_cn_with_ab_purity.txt")
inputs = c(inputs, "methylation_vs_cn_with_ab_purity_pam50.txt")

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

for(i in 1:4){
  
  inputi  = inputs[i]
  outputi = outputs[i]

  # ------------------------------------------------------------------------
  # read in eQTL results
  # ------------------------------------------------------------------------

  pvs = read.table(inputi, header=TRUE, sep="\t", as.is=TRUE)
  dim(pvs)
  pvs[1:2,]

  summary(pvs$p.value)
  summary(pvs$FDR)

  table(pvs$p.value < 1e-30)
  table(pvs$p.value < 1e-40)
  table(pvs$p.value < 1e-50)
  table(pvs$p.value < 1e-60)
  table(pvs$p.value < 1e-80)

  # ------------------------------------------------------------------------
  # plot it
  # ------------------------------------------------------------------------

  table(pvs$gene %in% infoE$Composite.Element.REF)
  table(pvs$SNP  %in% infoM$gene)

  ptplot   = pvs
  geneID   = match(ptplot$gene, infoE$Composite.Element.REF)
  markerID = match(ptplot$SNP,  infoM$gene)
  scores   = ptplot$p.value
  scuts    = c(1e-20, 1e-30, 1e-40, 1e-60)
  cols     = c("green", "blue", "red", "black")
  eChr     = infoE$Chromosome
  ePos     = infoE$Genomic_Coordinate
  mChr     = gsub("chr", "", infoM$chr)
  mPos     = 0.5*(infoM$start + infoM$end)
  chroms   = 1:22

  png(sprintf("../figures2/%s", outputi), width=7.5, height=9,
  res=200, units="in")
  
  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab="CN Location", ylab="Methylation Probe Location",
  plot.hotspots=TRUE, hotspots.cut=10, score.type="p-value")
  dev.off()

  # ------------------------------------------------------------
  # plot it after removing those possible problematic probes
  # ------------------------------------------------------------

  table(pvs$gene  %in% infMd$Name)
  ptplot   = pvs[which(pvs$gene  %in% infMd$Name),]

  geneID   = match(ptplot$gene, infoE$Composite.Element.REF)
  markerID = match(ptplot$SNP,  infoM$gene)
  scores   = ptplot$p.value
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




q(save = "no")

