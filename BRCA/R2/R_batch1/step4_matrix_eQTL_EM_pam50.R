
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# prepare new covariate file
# ------------------------------------------------------------

cDat = read.table(file="cov_EM_with_ab_purity_pam50.txt", sep = "\t",
                  header = TRUE, as.is=TRUE)

dim(cDat)
cDat[c(1:2,(nrow(cDat)-5):nrow(cDat)),1:5]

cDat = cDat[-which(cDat$id=="purity"),]

dim(cDat)
cDat[c(1:2,(nrow(cDat)-5):nrow(cDat)),1:5]

write.table(cDat, file = "cov_EM_with_pam50.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# prepare new covariate file
# ------------------------------------------------------------

library(MatrixEQTL)

SNP_file_name        = "methylation_mvalue.txt"
covariates_file_name = "cov_EM_with_pam50.txt"
expression_file_name = "expression_log_TReC.txt"
output_file_name     = "expression_vs_methylation_pam50.txt"

useModel              = modelLINEAR;
pvOutputThreshold     = 1e-40;

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

