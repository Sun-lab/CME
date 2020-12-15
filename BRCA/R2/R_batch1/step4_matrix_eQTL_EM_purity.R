
library(MatrixEQTL)

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

ff0 = "patient_brca_female_Caucasian_EMC_info_absolute.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

cDat = read.table(file="cov_EM.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(cDat)
cDat[c(1:2,(nrow(cDat)-2):nrow(cDat)),1:5]

table(names(cDat)[-1] == emInfo$patient_id)
purity = data.frame(id="purity", t(emInfo$abs_purity))
names(purity) = names(cDat)

cDat = rbind(cDat, purity)
dim(cDat)
cDat[c(1:3,(nrow(cDat)-4):nrow(cDat)),1:5]

write.table(cDat, file = "cov_EM_with_absolute_purity.txt",
append = FALSE, quote = FALSE, sep = "\t",
row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# run matrix eQTL
# ------------------------------------------------------------

SNP_file_name        = "methylation_mvalue.txt"
covariates_file_name = "cov_EM_with_absolute_purity.txt"
expression_file_name = "expression_log_TReC.txt"
output_file_name     = "expression_vs_methylation_absolute_purity.txt"

useModel              = modelLINEAR;
pvOutputThreshold     = 1e-30;

#errorCovariance = numeric();

snps = SlicedData$new(); 
snps$fileDelimiter = '\t'; # the TAB character 
snps$fileOmitCharacters = 'NA'; # denote missing values; 
snps$fileSkipRows = 1; # one row of column labels 
snps$fileSkipColumns = 1; # one column of row labels 
snps$fileSliceSize = 2000; # read file in pieces of 2000 rows 
snps$LoadFile( SNP_file_name );

genes = SlicedData$new(); 
genes$fileDelimiter = '\t'; # the TAB character 
genes$fileOmitCharacters = 'NA'; # denote missing values; 
genes$fileSkipRows = 1; # one row of column labels 
genes$fileSkipColumns = 1; # one column of row labels 
genes$fileSliceSize = 2000; # read file in pieces of 2000 rows 
genes$LoadFile( expression_file_name );


cvrt = SlicedData$new(); 
cvrt$fileDelimiter = '\t'; # the TAB character 
cvrt$fileOmitCharacters = 'NA'; # denote missing values; 
cvrt$fileSkipRows = 1; # one row of column labels 
cvrt$fileSkipColumns = 1; # one column of row labels 
cvrt$fileSliceSize = 2000; # read file in pieces of 2000 rows 
cvrt$LoadFile( covariates_file_name );

gc()

me = Matrix_eQTL_main(
snps, 
genes, 
cvrt, 
output_file_name = output_file_name, 
pvOutputThreshold = pvOutputThreshold,
useModel = modelLINEAR, 
errorCovariance = numeric(), 
verbose = TRUE, 
pvalue.hist = T
)

q(save = "no")

