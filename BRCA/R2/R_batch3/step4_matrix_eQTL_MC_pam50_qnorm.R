
library(MatrixEQTL)

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

SNP_file_name        = "cn_values.txt"
covariates_file_name = "cov_EM_with_pam50.txt"
expression_file_name = "methylation_mvalue_qnorm.txt"
output_file_name     = "methylation_vs_cn_with_pam50_qnorm.txt"

useModel              = modelLINEAR;
pvOutputThreshold     = 1e-20;

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

