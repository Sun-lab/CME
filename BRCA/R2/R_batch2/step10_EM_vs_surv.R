
library(survival)

# ------------------------------------------------------------
# read in location annotation
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA/_data2")

info = read.table("methylation_eQTM_info.txt", sep="\t", header=TRUE,
                    as.is=TRUE, quote="")
dim(info)
info[1:2,]

pval = read.table("methylation_eQTM_pval.txt", sep="\t", header=TRUE,
                    as.is=TRUE, quote="")

dim(pval)
pval[1:2,]

table(pval$SNP %in% info$Name)
pval = pval[which(pval$SNP %in% info$Name),]
dim(pval)

table(info$Genome_Build)

# ------------------------------------------------------------
# read in clinical information, check survival time after
# conditioning on stage and plate will make the model
# unidentifiable. so we will not consider these covariates.
# ------------------------------------------------------------

ff1  = "patient_brca_female_Caucasian_EMC_info_absolute.txt"

pInfo = read.table(ff1, sep = "\t", header=TRUE, as.is=TRUE)
dim(pInfo)
pInfo[1:2,]

table(pInfo$vital_status)
table(pInfo$vital_status, pInfo$death_days_to=="[Not Applicable]")

status = rep(0, nrow(pInfo))
status[which(pInfo$vital_status == "Dead")] = 1

time   = pInfo$last_contact_days_to
time[which(status==1)] = as.numeric(pInfo$death_days_to[which(status==1)])
summary(time)
time[time < 0] = 1

c1 = coxph(Surv(time, status) ~ pInfo$age_at_diagnosis)
c1

t1 = table(pInfo$ajcc_pathologic_tumor_stage, useNA="ifany")
t1
names(t1)

stage = rep(NA, nrow(pInfo))
stage[which(pInfo$ajcc_pathologic_tumor_stage %in% names(t1)[1:2])] = 1
stage[which(pInfo$ajcc_pathologic_tumor_stage %in% names(t1)[3:5])] = 2
stage[which(pInfo$ajcc_pathologic_tumor_stage %in% names(t1)[6:9])] = 3
stage[which(pInfo$ajcc_pathologic_tumor_stage %in% names(t1)[10])]  = 4

c1  = coxph(Surv(time, status) ~ stage)
c1

c1  = coxph(Surv(time, status) ~ tissue_source_site + pam50 + stage,
            data=pInfo)
c1

# ------------------------------------------------------------
# read in gene expression data and methylation data
# ------------------------------------------------------------

datE = read.table(file = "expression_log_TReC_rm_cn.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datE)
datE[1:2,1:5]

datM = read.table(file = "methylation_mvalue.txt", sep = "\t",
header = TRUE, as.is=TRUE)
dim(datM)
datM[1:2,1:5]

table(names(datE) == names(datM))
table(names(datE)[-1] == pInfo$patient_id)

# ------------------------------------------------------------
# evaluate the association between gene expression,
# DNA methylation and survival time
# ------------------------------------------------------------

pDatE  = data.matrix(datE[,-1])
survPE = rep(NA, nrow(datE))

for(i in 1:nrow(datE)){
  if(i %% 1000 == 0){
    cat(i, date(), "\n")
  }

  e1 = pDatE[i,]
  c1 = coxph(Surv(time, status) ~ e1)
  s1 = summary(c1)
  survPE[i] = s1$coef[5]
}


pDatM  = data.matrix(datM[,-1])
survPM = rep(NA, nrow(datM))

for(i in 1:nrow(datM)){
  
  if(i %% 10000 == 0){
    cat(i, date(), "\n")
  }
  
  m1 = pDatM[i,]
  c1 = coxph(Surv(time, status) ~ m1)
  s1 = summary(c1)
  survPM[i] = s1$coef[5]

}

# ------------------------------------------------------------
# functions of qq plots
# ------------------------------------------------------------


qqp <- function(pvals, main, confidence=.95, cutoff=1){
  
  alpha = 1-confidence
  n     = length(pvals)
  
  pvals[is.na(pvals)]=1
  pvals=sort(pvals)
  
  k=c(1:n)
  
  lower = cutoff*qbeta(alpha/2, k, n+1-k)
  upper = cutoff*qbeta((1-alpha/2), k, n+1-k)
  
  expected = cutoff*k/(n+1)
  
  biggest= max(-log10(pvals), -log10(expected))
  
  plot(-log10(expected), -log10(pvals), xlim=c(0,biggest),
  ylim=c(0,biggest), pch=20, xlab="-log10(expected p-value)",
  ylab="-log10(observed p-value)", cex=0.6, bty="n",
  main=main, col="darkgrey")
  
  lines(-log10(expected), -log10(lower), lty=2)
  lines(-log10(expected), -log10(upper), lty=2)
  
}

qqp.add <- function(pvals, cutoff=1, col="red"){
  
  n = length(pvals)
  
  pvals[is.na(pvals)]=1
  pvals=sort(pvals)
  
  k=c(1:n)
  
  expected = cutoff*k/(n+1)
  biggest= max(-log10(pvals), -log10(expected))
  
  points(-log10(expected), -log10(pvals), col=col, pch=9, cex=0.6)
  
}

# ------------------------------------------------------------
# check p-values of those eQTMs
# ------------------------------------------------------------

meths = unique(pval$SNP)
genes = unique(pval$gene)


surPE1 = survPE[match(genes, datE$id)]
surPM1 = survPM[match(meths, datM$id)]

length(surPE1)
length(surPM1)

png("../figures2/surv_eQTM.png", width=8, height=4, units="in", res=400)
par(mar=c(5,4,2,1), bty="n", mfrow=c(1,2))
qqp(survPE, "survival time vs. gene expression")
qqp.add(surPE1)
legend("bottomright", c("all genes", "genes of eQTMs"), bty="n",
  col=c("darkgrey", "red"), pch=c(20, 9))

qqp(survPM, "survival time vs. DNA methylation")
qqp.add(surPM1)
legend("bottomright", c("all methy probes", "eQTMs"), bty="n",
  col=c("darkgrey", "red"), pch=c(20, 9))
dev.off()

# ------------------------------------------------------------
# check p-values of those hot methy-probes
# ------------------------------------------------------------

mHot = read.table("methylation_hot.txt", sep="\t", header=T, as.is=T, quote="")
dim(mHot)
mHot[1:2,]

surPM1 = survPM[match(mHot$methyProbe, datM$id)]

gene1e60 = read.table("gene1e60.txt", header=TRUE, sep="\t", as.is=TRUE)
dim(gene1e60)
gene1e60[c(1:5,96:100,360:362),]

gene100 = gene1e60[gene1e60$freqIn1e60 >= 100,]
dim(gene100)
gene100[1:2,]

surPE1 = survPE[match(gene100$gene, datE$id)]


png("../figures2/surv_hotEM.png", width=8, height=4, units="in", res=400)
par(mar=c(5,4,2,1), bty="n", mfrow=c(1,2))
qqp(survPE, "survival time vs. gene expression")
qqp.add(surPE1)
legend("bottomright", c("all genes", "hot genes"), bty="n",
col=c("darkgrey", "red"), pch=c(20, 9))

qqp(survPM, "survival time vs. DNA methylation")
qqp.add(surPM1)
legend("bottomright", c("all methy probes", "hot methy probes"), bty="n",
col=c("darkgrey", "red"), pch=c(20, 9))
dev.off()


q(save = "no")

