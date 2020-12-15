
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

inputs = "expression_vs_methylation.txt"
inputs = c(inputs, "expression_vs_methylation_pam50.txt")
inputs = c(inputs, "expression_vs_methylation_absolute_purity.txt")
inputs = c(inputs, "expression_vs_methylation_ab_purity_pam50.txt")

inputs

outputs = gsub(".txt", ".png", inputs, fixed=TRUE)

outputs

# --------------------------------------------------------------
# read in location annotation
# --------------------------------------------------------------

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
  
  # ------------------------------------------------------------
  # read in eQTL results
  # ------------------------------------------------------------

  pvs = read.table(inputi, header=TRUE, sep="\t", as.is=TRUE)
  dim(pvs)
  pvs[1:2,]

  summary(pvs$p.value)

  pvs = pvs[pvs$p.value < 1e-50,]
  dim(pvs)
  gc()

  summary(pvs$p.value)

  table(pvs$p.value < 1e-60)
  table(pvs$p.value < 1e-70)
  table(pvs$p.value < 1e-80)
  table(pvs$p.value < 1e-90)
  table(pvs$p.value < 1e-100)

  # ------------------------------------------------------------
  # plot it
  # ------------------------------------------------------------

  table(pvs$gene %in% infoE$gene)
  table(pvs$SNP  %in% infoM$Composite.Element.REF)

  geneID   = match(pvs$gene, infoE$gene)
  markerID = match(pvs$SNP,  infoM$Composite.Element.REF)
  scores   = pvs$p.value
  scuts    = c(1e-60, 1e-70, 1e-80, 1e-100)
  cols     = c("green", "blue", "red", "black")
  eChr     = gsub("chr", "", infoE$chr)
  ePos     = 0.5*(infoE$start + infoE$end)
  mChr     = infoM$Chromosome
  mPos     = infoM$Genomic_Coordinate
  chroms   = 1:22

  png(sprintf("../figures2/%s", outputi), width=7.5, height=9,
      res=200, units="in")

  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
            mPos, chroms, xlab="Methylation Probe Location",
            ylab="Transcript Location", plot.hotspots=TRUE,
            hotspots.cut=10, score.type="p-value")
            
  dev.off()

  # ------------------------------------------------------------
  # plot it after removing those possible problematic probes
  # ------------------------------------------------------------

  table(pvs$SNP  %in% infMd$Name)
  pvs   = pvs[which(pvs$SNP  %in% infMd$Name),]

  geneID   = match(pvs$gene, infoE$gene)
  markerID = match(pvs$SNP,  infoM$Composite.Element.REF)
  scores   = pvs$p.value
  scuts    = c(1e-60, 1e-70, 1e-80, 1e-100)
  cols     = c("green", "blue", "red", "black")
  eChr     = gsub("chr", "", infoE$chr)
  ePos     = 0.5*(infoE$start + infoE$end)
  mChr     = infoM$Chromosome
  mPos     = infoM$Genomic_Coordinate
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

# ------------------------------------------------------------------------
# check the association between purity and subtypes
# ------------------------------------------------------------------------

ff0 = "patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

lm1 = lm(emInfo$abs_purity ~ emInfo$pam50)
summary(lm1)
anova(lm1)

pdf("../figures2/purity_vs_subtypes.pdf", width=5, height=3)
par(mar=c(2,4,1,1), bty="n")
boxplot(emInfo$abs_purity ~ emInfo$pam50)

stripchart(emInfo$abs_purity ~ emInfo$pam50,
method = "jitter", jitter = 0.2, offset = 1/3,
vertical = TRUE, add=TRUE, pch=1)
dev.off()

q(save = "no")

