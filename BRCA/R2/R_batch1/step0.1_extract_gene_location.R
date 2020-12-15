
# ------------------------------------------------------------
# read in the information of features
# ------------------------------------------------------------

setwd("~/research/TCGA/_Sun_MethyE/BRCA")

ff0   = "METADATA/UNC__IlluminaHiSeq_RNASeqV2/TCGA.hg19.June2011.gaf"
infoE = read.table(ff0, sep = "\t", header = TRUE, as.is=TRUE, comment.char="")
dim(infoE)
infoE[1:2,]
length(unique(infoE$FeatureID))

table(infoE$FeatureType, useNA="ifany")

infoE = infoE[infoE$FeatureType=="gene",]
dim(infoE)
infoE[1:2,]
length(unique(infoE$FeatureID))

table(infoE$Composite)
table(infoE$FeatureID == infoE$Gene)

cols = c("FeatureID", "Gene", "GeneLocus")
infoE[which(infoE$FeatureID != infoE$Gene), cols][1:5,]

# here I remove those genes that appear more than once...

# ------------------------------------------------------------
# take a subset of information for gene location
# ------------------------------------------------------------

infoE1 = infoE[which(infoE$FeatureID == infoE$Gene),16:19]
dim(infoE1)
infoE1[1:5,]

length(unique(infoE1$Gene))
infoE1 = unique(infoE1)
dim(infoE1)
infoE1[1:5,]

table(infoE1$FeatureAliases=="")
table(infoE1$FeatureInfo=="")

infoE1 = infoE1[,1:2]
dim(infoE1)
infoE1[1:5,]

gene1 = strsplit(infoE1$Gene, split="|", fixed=TRUE)
table(sapply(gene1, length))

gene1 = matrix(unlist(gene1), byrow=TRUE, ncol=2)
dim(gene1)
gene1[1:5,]

locus1 = strsplit(infoE1$GeneLocus, split=":", fixed=TRUE)
table(sapply(locus1, length))

locus1 = matrix(unlist(locus1), byrow=TRUE, ncol=3)
dim(locus1)
locus1[1:5,]

locus2 = strsplit(locus1[,2], split="-", fixed=TRUE)
table(sapply(locus2, length))

locus2 = matrix(unlist(locus2), byrow=TRUE, ncol=2)
dim(locus2)
locus2[1:5,]

infoE2 = data.frame(cbind(infoE1$Gene, gene1, locus1[,1], locus2, locus1[,3]),
                    stringsAsFactors=FALSE)
dim(infoE2)
infoE2[1:5,]

names(infoE2) = c("gene", "geneSymbol", "geneID", "chr", "start", "end", "strand")
dim(infoE2)
infoE2[1:5,]

table(infoE2$chr)
table(infoE2$strand)

infoE2 = infoE2[infoE2$chr %in% paste("chr", c(1:22,"X","Y"), sep=""),]
dim(infoE2)
infoE2[1:5,]

# ------------------------------------------------------------
# mapping from ensembl ID to gene ID
# ------------------------------------------------------------

g2e = read.table("~/research/data/human/hg19_info/gene2ensembl",
                  sep="\t", header=TRUE, as.is=TRUE)
dim(g2e)
g2e[1:2,]

g2e = g2e[g2e$tax_id==9606,]
dim(g2e)
g2e[1:2,]

table(infoE2$geneID %in% g2e$GeneID)

infoE2$ensembl = rep(NA, nrow(infoE2))

mat1 = match(infoE2$geneID, g2e$GeneID)
wmat = which(!is.na(mat1))

infoE2$ensembl[wmat] = g2e$Ensembl_gene_identifier[mat1[wmat]]

dim(infoE2)
infoE2[1:5,]

write.table(infoE2, file = "_data2/gene_info.txt", append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

q(save = "no")

