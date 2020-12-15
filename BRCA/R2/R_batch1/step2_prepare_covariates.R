
# ------------------------------------------------------------
# read in the information of the samples to be used
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA")

ff0    = "_data2/patient_brca_female_Caucasian_EMC_info.txt"
emInfo = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(emInfo)
emInfo[1,]

table(is.na(emInfo$cn_file))
colSums(is.na(emInfo))

emInfo = emInfo[which(!is.na(emInfo$cn_file)),]
dim(emInfo)
emInfo[1,]
colSums(is.na(emInfo))

# ------------------------------------------------------------
# read in purity information from pan12 study
# ------------------------------------------------------------

ff0    = "../shared_data/pancan12/pancan12.sample_info_filtered.txt"
purity = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE)
dim(purity)
purity[1:2,]

table(emInfo$cn_barcode %in% purity$tcga_id)
emInfo = emInfo[which(emInfo$cn_barcode %in% purity$tcga_id),]
dim(emInfo)
table(emInfo$cn_barcode %in% purity$tcga_id)

purity = purity[match(emInfo$cn_barcode, purity$tcga_id),]
dim(purity)
purity[1:2,]

table(emInfo$cn_barcode == purity$tcga_id)
table(purity$abs_call)

emInfo = cbind(emInfo, purity[,7:10])

# ------------------------------------------------------------
# read in subtype information
# ------------------------------------------------------------

ff1     = "_data2/ABSOLUTE_purity_and_ploidy_BRCA.txt"
subtype = read.table(ff1, sep = "\t", header = TRUE, as.is=TRUE)
dim(subtype)
subtype[1:2,]

table(emInfo$bcr_patient_barcode %in% subtype$individual_id)

subtype = subtype[match(emInfo$bcr_patient_barcode, subtype$individual_id),]
dim(subtype)
subtype[1:2,]
table(emInfo$bcr_patient_barcode == subtype$individual_id)

fig0  = "../shared_data/pancan12/pan12_purity_ploidy_checking.pdf"

pdf(fig0, width=8, height=4)
par(mar=c(5,4,1,1), mfrow=c(1,2), bty="n")
plot(emInfo$abs_purity, subtype$absolute_extract_purity)
plot(emInfo$abs_ploidy, subtype$absolute_extract_ploidy)
dev.off()

emInfo = cbind(emInfo, subtype$PAM50.Call)
dim(emInfo)
names(emInfo)[21] = "pam50"
emInfo[1:2,]

# ------------------------------------------------------------
# here we only keep those samples with purity and subtypes
# for the ease of comparision when we do eQTL mapping
# with or without purity estimates
# ------------------------------------------------------------

colSums(is.na(emInfo))

emInfo = emInfo[which(!is.na(emInfo$pam50)),]

dim(emInfo)
emInfo[1:2,]

colSums(is.na(emInfo))

# ------------------------------------------------------------
# obtain sample information of methylation data
# because the column of 'methylation_barcode' may contain the
# barcode for more than one sample, I used the column of
# methylation_file here.
# ------------------------------------------------------------

samM = strsplit(emInfo$methylation_file, split="-", fixed=TRUE)
table(sapply(samM, length))
samM = matrix(unlist(samM), byrow=TRUE, ncol=9)
dim(samM)
samM[1:2,]
table(samM[,1])
table(samM[,3])
table(samM[,9])

table(samM[,2], samM[,8])

samM = samM[,4:8]
dim(samM)
samM[1:2,]

samM = data.frame(samM, stringsAsFactors=FALSE)
names(samM) = c("institution", "patientID", "type", "portion", "plate")
dim(samM)
samM[1:2,]

length(unique(samM$patientID))
table(samM$institution == emInfo$tissue_source_site)

apply(samM[,-2], 2, function(v){sort(table(v))})

w2rm   = which(samM$institution=="GI" | samM$institution=="HN")
emInfo = emInfo[-w2rm,]

dim(emInfo)
emInfo[1:2,]
colSums(is.na(emInfo))

samM = samM[-w2rm,]
dim(samM)
samM[1:2,]

apply(samM[,-2], 2, function(v){sort(table(v))})

# ------------------------------------------------------------
# obtain sample information of gene expression data
# ------------------------------------------------------------

samE = strsplit(emInfo$expression_barcode, split="-", fixed=TRUE)
table(sapply(samE, length))
samE = matrix(unlist(samE), byrow=TRUE, ncol=7)
dim(samE)
samE[1,]

table(samE[,1])
table(samE[,7])

samE = data.frame(samE[,2:6], stringsAsFactors=FALSE)
names(samE) = c("institution", "patientID", "type", "portion", "plate")
dim(samE)
samE[1:2,]

apply(samE[,-2], 2, function(v){sort(table(v))})

table(samM$institution == samE$institution)
table(samM$patientID == samE$patientID)
table(samM$plate == samE$plate)

## plate ID of DNA methylation almost perfectly nested
## within the plateID of gene expression data, so we can
## ignore the DNA methylation plate IDs

table(samM$plate, samE$plate)
length(unique(samM$plate))
length(unique(samE$plate))

table(samM$institution)
table(samM$plate)

# ------------------------------------------------------------
# read in genotype PC information
# ------------------------------------------------------------

raceD = read.table("QC_PCA/final_data2_CaucasianOnly_BRCA.txt",
                  header=TRUE, sep="\t")
dim(raceD)
raceD[1:2,]

pid = substring(raceD$TCGA_ID, 1, 12)
table(emInfo$bcr_patient_barcode %in% pid)
raceD = raceD[match(emInfo$bcr_patient_barcode, pid),]
dim(raceD)
raceD[1:2,]

emInfo = cbind(emInfo, raceD[,3:5])
dim(emInfo)
emInfo[1:2,]

# ------------------------------------------------------------
# write out results
# ------------------------------------------------------------

ff1  = "_data2/patient_brca_female_Caucasian_EMC_info_absolute.txt"

write.table(emInfo, file = ff1, append = FALSE, quote = FALSE, sep = "\t",
row.names = FALSE, col.names = TRUE)


# ------------------------------------------------------------
# create the design matrix
# ------------------------------------------------------------

table(samM$patientID == emInfo$patient_id)

institute = as.factor(samE$institution)
plate     = as.factor(samE$plate)
age       = emInfo$age_at_diagnosis

sort(table(institute))
sort(table(plate))
length(unique(institute))
length(unique(plate))

cDat = model.matrix(~ institute + plate + age + noHM_PC1 + noHM_PC2 + noHM_PC3, data=emInfo)
dim(cDat)

cDat[1:2,]
cDat = cDat[,-1]

cor1 = cor(cDat)
eg1  = eigen(cor1)
eg1$values

colnames(cDat) = gsub("institute", "", colnames(cDat))
colnames(cDat) = gsub("plate", "", colnames(cDat))

cDat = t(cDat)
colnames(cDat) = emInfo$patient_id

dim(cDat)
cDat[1:5,1:5]

cDat = data.frame(id=rownames(cDat), cDat)

dim(cDat)
cDat[1:5,1:5]

# ------------------------------------------------------------
# Write out the results
# ------------------------------------------------------------

write.table(cDat, file = "_data2/cov_EM.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


q(save = "no")
