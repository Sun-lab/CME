
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

inputs = sprintf("expression_vs_methylation_rm_cn_PCs%d_qnorm.txt", 1:7)
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

labels = c("PC1", "PC1-2", "PC1-3", "PC1-4", "PC1-5", "PC1-6", "PC1-7")

pcuts = c(0, 10^(seq(-70, -30, by=10)))
np    = length(pcuts) - 1

tbl1 = tbl2 = matrix(nrow=length(inputs), ncol=np)

for(i in 1:length(inputs)){
  
  inputi  = inputs[i]
  outputi = outputs[i]
  
  # ------------------------------------------------------------
  # read in eQTL results
  # ------------------------------------------------------------
  
  pvs = read.table(inputi, header=TRUE, sep="\t", as.is=TRUE)
  dim(pvs)
  pvs[1:2,]
  
  table(pvs$SNP  %in% infMd$Name)
  pvs = pvs[which(pvs$SNP  %in% infMd$Name),]

  summary(pvs$p.value)
  
  peChr = eChr[match(pvs$gene, infoE$gene)]
  pmChr = mChr[match(pvs$SNP,  infoM$Composite.Element.REF)]

  for(j in 1:np){
    
    ww1 = which(pvs$p.value < pcuts[j+1] & pvs$p.value >= pcuts[j])
    tbl1[i,j] = length(ww1)
    tbl2[i,j] = length(which(peChr[ww1] == pmChr[ww1]))/length(ww1)
  }
  
  # ------------------------------------------------------------
  # plot it
  # ------------------------------------------------------------
  
  
  geneID   = match(pvs$gene, infoE$gene)
  markerID = match(pvs$SNP,  infoM$Composite.Element.REF)
  scores   = pvs$p.value
  cols     = c("green", "blue", "red", "black")

  if(i == 1){
    
    scuts    = c(1e-60, 1e-70, 1e-80, 1e-100)
    
    png("../figures2/expression_vs_methylation_ECM_purity1_stringent_qnorm.png",
        width=7.5, height=9, res=400, units="in")
    
    eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
    mPos, chroms, xlab="Methylation Probe Location",
    ylab="Transcript Location", plot.hotspots=TRUE,
    hotspots.cut=10, score.type="p-value")
    
    dev.off()
    
  }
  
  # ------------------------------------------------------------
  # plot it with more liberal cutoff
  # ------------------------------------------------------------
  
  scuts    = c(1e-30, 1e-40, 1e-50, 1e-60)
  
  outputi1 = sub(".png", "_rm_bad_m_qnorm.png", outputi, fixed=TRUE)
  
  png(sprintf("../figures2/liberal_%s", outputi1), width=7.5, height=9,
      res=400, units="in")
  
  eqtl.plot(geneID, markerID, scores, scuts, cols, eChr, ePos, mChr,
            mPos, chroms, xlab="Methylation Probe Location",
            ylab="Transcript Location", plot.hotspots=TRUE,
            hotspots.cut=10, score.type="p-value")
  
  dev.off()

}

pp1s = pcuts[1:(length(pcuts)-1)]
pp2s = pcuts[2:length(pcuts)]

rownames(tbl1) = rownames(tbl2) = labels
colnames(tbl1) = colnames(tbl2) = sprintf("[%.0e, %.0e)", pp1s, pp2s)

tbl1 = cbind(tbl1, rowSums(tbl1))
colnames(tbl1)[ncol(tbl1)] = "Total"

localPercent = rowSums(tbl1[,1:5]*tbl2[,1:5])/tbl1[,6]
tbl2 = cbind(tbl2, localPercent)

tbl2 = round(tbl2, 3)

tbl1
tbl2

tbl3 = matrix(nrow=7, ncol=6)

for(i in 1:nrow(tbl3)){
  for(j in 1:ncol(tbl3)){
    tbl3[i,j] = sprintf("%d (%.2f)", tbl1[i,j], tbl2[i,j])
  }
}

rownames(tbl3) = rownames(tbl1)
colnames(tbl3) = colnames(tbl1)
write.table(tbl3, sep=" & ", eol = "\\\\\n", quote=FALSE)

write.table(tbl1, file = "eQTL_with_ECM-purity_hits_qnorm.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = TRUE,
col.names = TRUE)

write.table(tbl2, file = "eQTL_with_ECM-purity_local_percentage_qnorm.txt",
  append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE,
col.names = TRUE)


q(save = "no")

