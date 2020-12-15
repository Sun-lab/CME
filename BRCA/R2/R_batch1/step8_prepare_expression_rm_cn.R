
# ------------------------------------------------------------
# read in gene expression and copy number data
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

eDat = read.table("expression_log_TReC.txt", sep = "\t",
                  header = TRUE, as.is=TRUE)
dim(eDat)
eDat[1:2,1:5]

cDat = read.table("cn_values.txt", sep = "\t", header = TRUE, as.is=TRUE)

dim(cDat)
cDat[1:2,1:5]

geneIds = intersect(eDat$id, cDat$id)
length(geneIds)

cDat = cDat[match(geneIds, cDat$id),]
eDat = eDat[match(geneIds, eDat$id),]

dim(eDat)
eDat[1:2,1:5]

dim(cDat)
cDat[1:2,1:5]

# ------------------------------------------------------------------------
# take residuals
# ------------------------------------------------------------------------

ePdat = data.matrix(eDat[,-1])
cPdat = data.matrix(cDat[,-1])

dim(ePdat)
ePdat[1:2,1:5]

dim(cPdat)
cPdat[1:2,1:5]

eNA = rowSums(is.na(ePdat))
cNA = rowSums(is.na(cPdat))

table(eNA)
table(cNA)

rDat = matrix(NA, nrow=nrow(ePdat), ncol=ncol(ePdat))
dim(rDat)

for(i in 1:nrow(ePdat)){
  if(i %% 2000 == 0){ cat(i, date(), "\n") }
  yi = ePdat[i,]
  xi = cPdat[i,]
  
  ri = lm(yi ~ xi)$resid
  rDat[i,which(!is.na(xi))] = ri
}

rNA = rowSums(is.na(rDat))
table(rNA)

rDatDF = data.frame(id=eDat$id, rDat)
names(rDatDF) = names(eDat)
dim(rDatDF)
rDatDF[1:2,1:5]


write.table(rDatDF, file = "expression_log_TReC_rm_cn.txt",
            append = FALSE, quote = FALSE, sep = "\t",
            row.names = FALSE, col.names = TRUE)


q(save = "no")

