
normscore = function(vec) {
  len  = length(na.omit(vec))+1
  rank = rank(na.omit(vec))
  ties = (rank - floor(rank)) > 0
  new.vec = vec[!is.na(vec)]
  new.vec[!ties]=qnorm(rank[!ties]/len)
  new.vec[ties] =0.5*(qnorm((rank[ties]+0.5)/len)+qnorm((rank[ties]-0.5)/len))
  vec[!is.na(vec)] = new.vec
  vec
}

# ----------------------------------------------------------------------
# read in data
# ----------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

datE = read.table(file = "expression_log_TReC.txt",
sep = "\t", header = TRUE, as.is = TRUE)
dim(datE)
datE[1:2, 1:5]


cDat = read.table(file="cov_EM.txt", sep = "\t",
header = TRUE, as.is = TRUE)

dim(cDat)
cDat[1:5, 1:5]
cDat$id

X = data.matrix(t(cDat[, -1]))
dim(X)
X[1:5, 1:5]

# ----------------------------------------------------------------------
# find wheher some genes that have outliers
# ----------------------------------------------------------------------

pDat = data.matrix(datE[,-1])
dim(pDat)


# ----------------------------------------------------------------------
# (95% - median)/IQR
# ----------------------------------------------------------------------

probs = c(0, 0.25, 0.5, 0.75, 1)
qs = apply(pDat, 1, quantile, probs = probs, na.rm=TRUE)

dim(qs)
qs = t(qs)
qs[1:5,]

summary(qs)

IQR = qs[,4] - qs[,2]
qs2 = qs[,c(1,5)]

qs2 = (qs2 - qs[,3])/IQR
dim(qs2)
qs2[1:2,]

outlier = apply(abs(qs2), 1, max)

summary(outlier)

table(outlier > 5)
table(outlier > 6)
table(outlier > 7)
table(outlier > 8)
table(outlier > 10)
table(outlier > 15)

# ----------------------------------------------------------------------
# quantile transform expression data
# ----------------------------------------------------------------------

pDatN = t(apply(pDat, 1, normscore))
dim(pDatN)
pDatN[1:2,1:5]

summary(apply(pDatN, 1, mean, na.rm=T))
summary(apply(pDatN, 1, sd, na.rm=T))

pDatN = data.frame(id=datE$id, pDatN)
dim(pDatN)
pDatN[1:2,1:5]

write.table(pDatN, file = "expression_log_TReC_qnorm.txt", append = FALSE,
  quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


# ----------------------------------------------------------------------
# normalize data after removing copy number effects
# ----------------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

datE = read.table(file = "expression_log_TReC_rm_cn.txt",
sep = "\t", header = TRUE, as.is = TRUE)
dim(datE)
datE[1:2, 1:5]

pDat = data.matrix(datE[,-1])
dim(pDat)

pDatN = t(apply(pDat, 1, normscore))
dim(pDatN)
pDatN[1:2,1:5]

summary(apply(pDatN, 1, mean, na.rm=T))
summary(apply(pDatN, 1, sd, na.rm=T))

pDatN = data.frame(id=datE$id, pDatN)
dim(pDatN)
pDatN[1:2,1:5]

write.table(pDatN, file = "expression_log_TReC_rm_cn_qnorm.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


q(save="no")
