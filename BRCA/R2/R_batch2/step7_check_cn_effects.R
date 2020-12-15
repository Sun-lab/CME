
setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

# ------------------------------------------------------------
# read in data
# ------------------------------------------------------------

cDat = read.table("cn_values.txt", header=TRUE, as.is=TRUE, sep="\t")
dim(cDat)
cDat[1:2,1:5]

infoC = read.table("cn_info.txt", sep="\t", header=TRUE, as.is=TRUE)
dim(infoC)
infoC[1:2,]

xDat = read.table("cov_EM_with_ab_purity_pam50.txt", header=TRUE,
                  as.is=TRUE, sep="\t")

dim(xDat)
xDat[1:2,1:5]

table(names(cDat) == names(xDat))

ids  = xDat$id
xDat = t(xDat[,-1])
xDat = data.frame(xDat)
names(xDat) = ids
dim(xDat)
xDat[1:2,1:5]

names(xDat)

X0 = data.matrix(xDat[,1:37])
X1 = data.matrix(xDat)

# ------------------------------------------------------------
# check assocaitoin between copy number and subtype
# ------------------------------------------------------------

cPdat = data.matrix(cDat[,-1])
dim(cPdat)

pvs = rep(NA, nrow(cPdat))

for(i in 1:length(pvs)){
  if(i %% 1000 == 0){
    cat(i, date(), "\n")
  }
  yi   = cPdat[i,]
  lmi0 = lm(yi ~ X0)
  lmi1 = lm(yi ~ X1)
  an1  = anova(lmi0, lmi1)
  
  pvs[i] = an1[2,6]
}

setwd("~/research/TCGA/_Sun_MethyE/BRCA/")

source("../shared_code/manhattan.R")

png("figures2/cn_vs_PAM50.png", width=10, height=3.5, res=400, units="in")
par(mar=c(5,4,1,1), bty="n")
manhattan(pvs, gsub("chr", "", infoC$chr), infoC$start)
dev.off()

q(save = "no")

