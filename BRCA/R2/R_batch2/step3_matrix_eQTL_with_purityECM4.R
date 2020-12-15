
library(MatrixEQTL)

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# prepare new covariate file
# ------------------------------------------------------------

cDat = read.table(file="cov_EM_with_ECM_purity3.txt", sep = "\t",
                  header = TRUE, as.is=TRUE)

dim(cDat)
cDat[c(1:2,(nrow(cDat)-2):nrow(cDat)),1:5]
cDat$id

pcs = read.table("ECM_pam50_eigen_vectors.txt", sep = "\t",
                  header = TRUE, as.is=TRUE)
dim(pcs)
pcs[1,]

cpDat = data.matrix(cDat[,-1])
cor(pcs[,1:3], t(cpDat[41:43,]))

pt1 = data.frame(id="purity_ECM_PC4", t(pcs[,4]))

names(pt1) = names(cDat)
dim(pt1)

cDat4 = rbind(cDat,  pt1)

dim(cDat4)
cDat4[(nrow(cDat4)-3):nrow(cDat4),1:5]
cDat4$id

write.table(cDat4, file = "cov_EM_with_ECM_purity4.txt",
append = FALSE, quote = FALSE, sep = "\t",
row.names = FALSE, col.names = TRUE)

# ------------------------------------------------------------
# run matrix eQTL
# ------------------------------------------------------------

SNP_file_name        = "methylation_mvalue.txt"
covariates_file_name = "cov_EM_with_ECM_purity4.txt"
expression_file_name = "expression_log_TReC_rm_cn.txt"
output_file_name     = "expression_vs_methylation_ECM_purity4.txt"

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

