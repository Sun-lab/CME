
# ------------------------------------------------------------
# read in the information of the samples to be used
# ------------------------------------------------------------

setwd("/lustre/scr/w/e/weisun/_Sun_methylE/BRCA/")

ff0    = "_data2/patient_brca_female_Caucasian_EM_info.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

# ------------------------------------------------------------
# read in methylation data
# ------------------------------------------------------------

setwd("./DNA_Methylation")

fs2use    = emInfo$methylation_file
ids2use   = emInfo$methylation_barcode

betaValue = matrix(NA, nrow=485577, ncol=length(fs2use))
colnames(betaValue) = ids2use

for(i in 1:length(fs2use)){
  
  if(i %% 10 ==0){
    cat(i, date(), "\n")
  }
  
  fi = paste("JHU_USC__HumanMethylation450/Level_3/", fs2use[i], sep="")
  d1 = scan(fi, what = character(), nlines=1)
  
  idsI = unlist(strsplit(ids2use[i], split="/", fixed=TRUE))
  
  if(any(! d1 %in% c("Hybridization", "REF", idsI))){
    stop("hi, something is wrong here.\n")
  }
  
  di = read.table(fi, sep="\t", header=TRUE, as.is=TRUE, skip=1)
  
  if(nrow(di) != 485577){
    stop("nrow(di) is not what we expected...\n")
  }
  
  if(i == 1){
    symbols = di$Composite.Element.REF
    rownames(betaValue) = symbols

    if(length(symbols) != length(unique(symbols))){
      stop("symbols are not unique\n")
    }
    info = di[,-2]
    dim(info)
    
    write.table(info, file = "HM450_info.txt", quote = FALSE,
      sep = "\t", row.names = FALSE, col.names = TRUE)
    
  }else{
    if(any(symbols != di$Composite.Element.REF)){
      stop("symbols do not match\n")
    }
  }
  
  betaValue[,i] = di$Beta_value
}

gc()

# ------------------------------------------------------------
# check the range of beta values,
# beta values are usually not extremely close to 0 or 1
# so we can calcuatle mvalues
# ------------------------------------------------------------

dim(betaValue)

nNAs = rowSums(is.na(betaValue))
w2kp = which(nNAs <= 0.2*ncol(betaValue))
betaValue = betaValue[w2kp,]
dim(betaValue)

gc()

sms = apply(betaValue, 2, summary)
dim(sms)
sms[,1:2]
apply(sms, 1, summary)

pdf("figures2/beta_value_1stQu_vs_median.pdf", width=5, height=5)
par(mar=c(5,4,1,1), bty="n")
plot(sms[2,], sms[3,], xlab="1st Qu.", ylab="Median")
dev.off()

colnames(sms)[sms[2,] > 0.1]

mValue = log(betaValue/(1-betaValue))
dim(mValue)

betaValue = signif(betaValue, 6)
mValue    = signif(mValue, 6)

write.table(betaValue, file = "betaValues_tumor_v2.txt",
  append = FALSE, quote = FALSE, sep = "\t",
  row.names = TRUE, col.names = TRUE)

write.table(mValue, file = "mValues_tumor_v2.txt",
  append = FALSE, quote = FALSE, sep = "\t",
  row.names = TRUE, col.names = TRUE)

rm(betaValue)
rm(mValue)

gc()

# ------------------------------------------------------------
# read in gene expression data
# ------------------------------------------------------------

setwd("../RNASeqV2")

fs2use  = emInfo$expression_file
ids2use = emInfo$expression_barcode

counts = matrix(NA, nrow=20531, ncol=length(fs2use))

for(i in 1:length(fs2use)){
  
  if(i %% 10 ==0){
    cat(i, date(), "\n")
  }
  
  fi = paste("UNC__IlluminaHiSeq_RNASeqV2/Level_3/", fs2use[i], sep="")
  di = read.table(fi, sep="\t", header=TRUE, as.is=TRUE)
  
  if(nrow(di) != 20531){
    stop("nrow(di) is not what we expected...\n")
  }

  if(i == 1){
    geneids = di$gene_id
  }else{
    if(any(di$gene_id != geneids)){
      stop("gene ids do not match\n")
    }
  }
  
  counts[,i] = di$raw_count
}

dim(counts)
counts[1:2,1:9]

rownames(counts) = geneids
colnames(counts) = ids2use

dim(counts)
counts[1:2,1:9]

write.table(counts, file = "rawCounts_tumor_v2.txt",
append = FALSE, quote = FALSE, sep = "\t",
row.names = TRUE, col.names = TRUE)

q(save = "no")

