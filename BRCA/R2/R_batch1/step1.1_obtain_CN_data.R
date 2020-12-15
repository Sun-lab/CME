
library(GenomicRanges)

# ------------------------------------------------------------
# read in the information of the samples to be used
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA")

ff0    = "_data2/patient_brca_female_Caucasian_EM_info.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

barcodes = strsplit(emInfo$expression_barcode, split="-")
table(sapply(barcodes, length))
barcodes = matrix(unlist(barcodes), ncol=7, byrow=TRUE)
dim(barcodes)
barcodes[1:2,]

barcodeE = apply(barcodes[,1:4], 1, paste, collapse="-")
barcodeE[1:5]

# ------------------------------------------------------------
# read in copy number file information
# ------------------------------------------------------------

ff1 = "METADATA/BI__Genome_Wide_SNP_6/broad.mit.edu_BRCA.Genome_Wide_SNP_6.sdrf.txt"
cnFile = read.table(ff1, sep = "\t", header = TRUE, as.is=TRUE)
dim(cnFile)
cnFile[1:5,1:5]
table(cnFile$Protocol.REF)

cnFile$Derived.Array.Data.File.3[1:5]
gl = grepl("nocnv_hg19.seg", cnFile$Derived.Array.Data.File.3)
table(gl)

table(cnFile$Derived.Array.Data.File.3[which(!gl)])
gl1 = grepl("hg19.seg", cnFile$Derived.Array.Data.File.1)
table(gl, gl1)

cnFile = cnFile[which(gl),]
cnFile = cnFile[,c("Comment..TCGA.Barcode.", "Derived.Array.Data.File.3")]
dim(cnFile)
cnFile[1:5,]

barcodes = strsplit(cnFile$Comment..TCGA.Barcode., split="-")
table(sapply(barcodes, length))
barcodes = matrix(unlist(barcodes), ncol=7, byrow=TRUE)
dim(barcodes)
barcodes[1:2,]

barcodeCN = apply(barcodes[,1:4], 1, paste, collapse="-")
barcodeCN[1:5]

table(barcodeE %in% barcodeCN)
mat1 = match(barcodeE, barcodeCN)

wNonNA = which(!is.na(mat1))
cnFileName = cnBarcode = rep(NA, nrow(emInfo))
cnFileName[wNonNA] = cnFile$Derived.Array.Data.File.3[mat1[wNonNA]]
cnBarcode[wNonNA]  = cnFile$Comment..TCGA.Barcode.[mat1[wNonNA]]

emInfo$cn_barcode = cnBarcode
emInfo$cn_file    = cnFileName

dim(emInfo)
emInfo[1,]

ff2 = "_data2/patient_brca_female_Caucasian_EMC_info.txt"

write.table(emInfo, file = ff2, append = FALSE, quote = FALSE, sep = "\t",
            eol = "\n", row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# read in gene location information
# ------------------------------------------------------------

ff3   = "_data2/gene_info.txt"
infoE = read.table(ff3, sep = "\t", header = TRUE, as.is=TRUE)
dim(infoE)
infoE[1:2,]
length(unique(infoE$gene))

table(infoE$chr, useNA="ifany")
table(infoE$strand, useNA="ifany")
summary(infoE$start)

gr1 = makeGRangesFromDataFrame(infoE, ignore.strand=TRUE)
gr1

# ------------------------------------------------------------
# read in copy number data
# ------------------------------------------------------------

fs2use    = emInfo$cn_file
ids2use   = emInfo$cn_barcode

cnValue = matrix(NA, nrow=nrow(infoE), ncol=nrow(emInfo))

for(i in 1:length(cnFileName)){
  
  if(i %% 10 ==0){
    cat(i, date(), "\n")
  }
  
  if(is.na(fs2use[i])){
    next
  }
  
  fi = paste("CNV_SNP_Array/BI__Genome_Wide_SNP_6/Level_3/", fs2use[i], sep="")
  
  if(! file.exists(fi)){
    next
  }
  
  d1 = read.table(fi, header=TRUE, as.is=TRUE, sep="\t")
  dim(d1)
  d1[1:2,]
  
  table(d1$Chromosome)
  d1$Chromosome = paste("chr", d1$Chromosome, sep="")
  
  gr2 = makeGRangesFromDataFrame(d1, ignore.strand=TRUE,
                                seqnames.field="Chromosome",
                                start.field="Start",
                                end.field="End")
  
  
  mtch = findOverlaps(gr1, gr2, select="first")
  table(is.na(mtch))
  
  ## the copy number information is often missing for genes at
  ## the begining or the end of the chromosome
  if(i == 1){
    cat("information of the genes without location information\n")
    eMiss = infoE[which(is.na(mtch)),]
    print(eMiss[order(eMiss$chr, eMiss$start),])
  }
  
  wNonNA = which(!is.na(mtch))
  cnValue[wNonNA,i] = d1$Segment_Mean[mtch[wNonNA]]
}

dim(cnValue)
cnValue[1:2,1:5]

colnames(cnValue) = emInfo$patient_id
rownames(cnValue) = infoE$gene
dim(cnValue)
cnValue[1:2,1:5]

table(rowSums(is.na(cnValue)))
0.2*ncol(cnValue)

w2kp = which(rowSums(is.na(cnValue)) < 0.2*ncol(cnValue))
cnValue = cnValue[w2kp,]

dim(cnValue)
cnValue[1:2,1:5]

table(colSums(is.na(cnValue)))
0.1*nrow(cnValue)

w2kp = which(colSums(is.na(cnValue)) < 0.1*nrow(cnValue))
cnValue = cnValue[,w2kp]

dim(cnValue)
cnValue[1:2,1:5]

write.table(cnValue, file = "CNV_SNP_Array/cn_data.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


q(save = "no")

