
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

inputs = "expression_vs_cn.txt"
inputs = c(inputs, "expression_vs_cn_with_pam50.txt")
inputs = c(inputs, "expression_vs_cn_with_absolute_purity.txt")
inputs = c(inputs, "expression_vs_cn_with_ab_purity_pam50.txt")

inputs

outputs = gsub(".txt", ".png", inputs, fixed=TRUE)

outputs

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

infoE = read.table("expression_info.txt", sep="\t", header=TRUE,
                  as.is=TRUE)
dim(infoE)
infoE[1:2,]

table(infoE$chr, useNA="ifany")

infoC = read.table("cn_info.txt", sep="\t", header=TRUE,
                  as.is=TRUE)
dim(infoC)
infoC[1:2,]
table(infoC$chr, useNA="ifany")

# --------------------------------------------------------------
# iteratively summarizing results
# --------------------------------------------------------------

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

  table(pvs$p.value < 1e-30)
  table(pvs$p.value < 1e-40)
  table(pvs$p.value < 1e-50)
  table(pvs$p.value < 1e-60)
  table(pvs$p.value < 1e-80)

  # ------------------------------------------------------------------------
  # plot it
  # ------------------------------------------------------------------------

  table(pvs$gene %in% infoE$gene)
  table(pvs$SNP  %in% infoC$gene)

  ptplot   = pvs
  geneID   = match(ptplot$gene, infoE$gene)
  markerID = match(ptplot$SNP,  infoC$gene)
  scores   = ptplot$p.value
  scuts    = c(1e-20, 1e-30, 1e-40, 1e-60)
  cols     = c("green", "blue", "red", "black")
  eChr     = gsub("chr", "", infoE$chr)
  ePos     = 0.5*(infoE$start + infoE$end)
  mChr     = gsub("chr", "", infoC$chr)
  mPos     = 0.5*(infoC$start + infoC$end)
  chroms   = 1:22

  png(sprintf("../figures2/%s", outputi), width=7.5, height=9,
  res=200, units="in")
  
  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
  mPos, chroms, xlab="CN Location", ylab="Transcript Location",
  plot.hotspots=TRUE, hotspots.cut=10, score.type="p-value")
  dev.off()

}


q(save = "no")

