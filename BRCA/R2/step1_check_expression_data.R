
library(data.table)
# ------------------------------------------------------------------------
# read in more detailed methyaltion information
# ------------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

datE = fread("expression_log_TReC.txt")
dim(datE)
datE[1:2,1:5]


datE2 = fread("expression_log_TReC_rm_cn_qnorm.txt")
dim(datE2)
datE2[1:2,1:5]

table(names(datE) == names(datE2))

length(intersect(datE$id, datE2$id))

samsE = read.table("patient_brca_female_Caucasian_EMC_info_absolute.txt", 
                   header=TRUE, sep="\t", as.is=TRUE)
dim(samsE)
samsE[1:2,]

table(names(datE)[-1] == samsE$patient_id)

q(save = "no")

