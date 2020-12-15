
# setwd("/home/groups/projects/ith/_Sun_methylE/BRCA/")

setwd("~/research/TCGA/_Sun_MethyE/BRCA/")

# ------------------------------------------------------------
# read in race information obtained by PCA
# ------------------------------------------------------------

raceD = read.table("QC_PCA/final_data_BRCA.txt", header=TRUE, sep="\t")
dim(raceD)
raceD[1:2,]

table(raceD$Caucasian, raceD$PASS.QC.1 + 2*raceD$PASS.QC.2, useNA="ifany")

pid = substring(raceD$TCGA_ID, 1, 12)

# ------------------------------------------------------------
# read in patient information
# ------------------------------------------------------------

ff0 = "Clinical/Biotab/nationwidechildrens.org_clinical_patient_brca.txt"
pat = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE,
                 na.string="[Not Available]", quote = "")

dim(pat)
pat[1:3,1:10]

# the first two rows are description, instead of data
pat = pat[-(1:2),]
dim(pat)
pat[1:2,1:10]
names(pat)

length(unique(pat$patient_id))

table(pat$prospective_collection, useNA="ifany")
table(pat$race, useNA="ifany")
table(pat$gender, useNA="ifany")
table(pat$ethnicity, useNA="ifany")

table(pid %in% pat$bcr_patient_barcode)
table(pat$bcr_patient_barcode %in% pid)

barcode = intersect(pid, pat$bcr_patient_barcode)
length(barcode)

raceD = raceD[match(barcode,pid),]
pat   = pat[match(barcode,pat$bcr_patient_barcode),]

dim(raceD)
dim(pat)
table(substring(raceD$TCGA_ID, 1, 12) == pat$bcr_patient_barcode)

table(pat$race, raceD$Caucasian, useNA="ifany")
table(pat$race, pat$ethnicity, useNA="ifany")
table(raceD$Caucasian, pat$ethnicity, useNA="ifany")

# ------------------------------------------------------------
# select thsoe female patients who are non-hispanic white
# ------------------------------------------------------------

pat1 = pat[which(pat$gender=="FEMALE" & raceD$Caucasian=="YES"), ]
dim(pat1)

table(pat1$tumor_status)
table(pat1$tumor_tissue_site)
table(pat1$tissue_source_site)
table(pat1$clinical_stage)
table(pat1$ajcc_pathologic_tumor_stage)

ff1 = "_data2/patient_brca_female_Caucasian.txt"
write.table(pat1, file=ff1, sep="\t", quote=FALSE, row.names = FALSE,
  col.names = TRUE)

nms  = c("bcr_patient_barcode", "patient_id", "tissue_source_site")
nms  = c(nms, "birth_days_to", "last_contact_days_to", "death_days_to")
nms  = c(nms, "vital_status", "ajcc_pathologic_tumor_stage", "tumor_status")
nms  = c(nms, "age_at_diagnosis")

pat2 = pat1[,nms]
dim(pat2)
pat2[1:2,]

ff2 = "_data2/patient_brca_female_Caucasian_short_table.txt"
write.table(pat2, file=ff2, sep="\t", quote=FALSE, row.names = FALSE,
  col.names = TRUE)

# ------------------------------------------------------------
# obtain the patient barcode for those with methylation data
# ------------------------------------------------------------

fls = read.table("file_manifest.txt", sep = "\t", header = TRUE, as.is=TRUE)

table(fls$Platform.Type, useNA="ifany")

flsDM = fls[fls$Platform.Type == "DNA Methylation",]
dim(flsDM)
flsDM[1:2,]

table(flsDM$Center)
table(flsDM$Platform)
table(flsDM$Level)

mids = strsplit(flsDM$Sample, split="-")
table(sapply(mids, length))
mids = matrix(unlist(mids), byrow=TRUE, ncol=4)
mids[1:2,]

mids = data.frame(mids, stringsAsFactors=TRUE)
names(mids) = c("tcga", "site", "patient", "sample")
table(mids$sample)

flsDM = flsDM[which(mids$sample=="01"),]

# ------------------------------------------------------------
# obtain the patient barcode for those with expression data
# ------------------------------------------------------------

flsEX = fls[which(fls$Platform.Type == "RNASeqV2"),]
dim(flsEX)
flsEX[1:2,]

flsEX = flsEX[grep("rsem.genes.results", flsEX$File.Name),]
dim(flsEX)
flsEX[1:2,]

table(flsEX$Center)
table(flsEX$Platform)
table(flsEX$Level)

table(table(flsEX$Barcode))

eids = strsplit(flsEX$Sample, split="-")
table(sapply(eids, length))
eids = matrix(unlist(eids), byrow=TRUE, ncol=4)
eids[1:2,]

eids = data.frame(eids, stringsAsFactors=TRUE)
names(eids) = c("tcga", "site", "patient", "sample")
table(eids$sample)

flsEX = flsEX[which(eids$sample=="01"),]

# ------------------------------------------------------------
# now we need to make sure the samples are from the same vial
# e.g, 01A matched with 01A instead of 01B
# ------------------------------------------------------------
# for methylation data, the column Barcode is not useable
# because it includes the barcode for more than one sample
# we will get barcode from file name
# ------------------------------------------------------------

barcodeDM = strsplit(flsDM$File.Name, split="-")
table(sapply(barcodeDM, length))

if(any(sapply(barcodeDM, length)!= 9)){
  stop("file name is not expected\n")
}

barcodeDM = matrix(unlist(barcodeDM), byrow=TRUE, ncol=9)
dim(barcodeDM)
barcodeDM[1:2,]

apply(barcodeDM, 2, table)

barcodeDM = cbind(rep("TCGA", nrow(barcodeDM)), barcodeDM[,4:6])
barcodeDM = apply(barcodeDM, 1, paste, collapse="-")
barcodeDM[1:5]

# ------------------------------------------------------------
# for gene expression data, we assume the sample are unique
# ------------------------------------------------------------

barcodeEX = strsplit(flsEX$Barcode, split="-")
table(sapply(barcodeEX, length))

if(any(sapply(barcodeEX, length)!= 7)){
  stop("barcodeEX is not expected\n")
}

barcodeEX = matrix(unlist(barcodeEX), byrow=TRUE, ncol=7)
dim(barcodeEX)
barcodeEX[1:2,]
barcodeEX = apply(barcodeEX[,1:4], 1, paste, collapse="-")
barcodeEX[1:2]

# ------------------------------------------------------------
# take intersection of expression and methyaltion data
# ------------------------------------------------------------

barcodes = sort(intersect(barcodeEX, barcodeDM))

length(barcodeEX)
length(barcodeDM)
length(barcodes)

flsDM = flsDM[match(barcodes, barcodeDM),]
flsEX = flsEX[match(barcodes, barcodeEX),]

dim(flsDM)
dim(flsEX)

flsDM[1:5,]
flsEX[1:5,]

if(any(flsDM$Sample != flsEX$Sample)){
  stop("sample name mismatch\n")
}

# ------------------------------------------------------------
# intersect with the samples with demographic/clinical data
# ------------------------------------------------------------

patientsEM = substring(barcodes,1,12)
patientsEM[1:2]
length(patientsEM)

patients2use = sort(intersect(patientsEM, pat2$bcr_patient_barcode))
length(patients2use)

flsDM2use = flsDM[match(patients2use, patientsEM),]
flsEX2use = flsEX[match(patients2use, patientsEM),]

dim(flsDM2use)
flsDM2use[1:2,]

dim(flsEX2use)
flsEX2use[1:2,]

pat3 = pat2[match(patients2use, pat2$bcr_patient_barcode),]
dim(pat3)
pat3[1:2,]

table(patients2use == pat3$bcr_patient_barcode)
table(patients2use == substr(flsDM2use$Sample, 1, 12))
table(patients2use == substr(flsEX2use$Sample, 1, 12))

pat3$methylation_barcode = flsDM2use$Barcode
pat3$methylation_file    = flsDM2use$File.Name

pat3$expression_barcode  = flsEX2use$Barcode
pat3$expression_file     = flsEX2use$File.Name

dim(pat3)
pat3[1:2,]

ff3 = "_data2/patient_brca_female_Caucasian_EM_info.txt"

write.table(pat3, file=ff3, sep="\t", quote=FALSE, row.names = FALSE,
col.names = TRUE)

q(save = "no")

