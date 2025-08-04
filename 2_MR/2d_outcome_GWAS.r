################
##Outcome GWAS##
################

# load libraries
library(SeqArray)
library(SeqArray)
library(SeqVarTools)
library(GENESIS)
library(SNPRelate)
library(Biobase)
library(ggplot2)
library(qqman)
library(GenomicRanges)
library(plyr); library(dplyr)


# convert VCF to GDS

vcffile <- "Genotyping.vcf.gz"
gdsfile <- "Genotyping.gds"

seqVCF2GDS(vcffile, gdsfile, fmt.import="GT", storage.option="LZMA_RA")


# load gds
gdsall <- seqOpen(gdsfile)

# run LD pruning
set.seed(100)
snpset <- snpgdsLDpruning(gdsall, method="corr", slide.max.bp=10e6, ld.threshold=sqrt(0.02),maf=0.01)

pruned <- unlist(snpset, use.names=FALSE)
length(pruned)


# compute the GRM
grm <- snpgdsGRM(gdsall, method="GCTA", snp.id = pruned)
names(grm)
dim(grm$grm)

# run KING-robust to get initial kinship estimates using KING-robust, 

king <- snpgdsIBDKING(gdsall, snp.id=pruned)
names(king)
dim(king$kinship)
kingMat <- king$kinship
colnames(kingMat) <- rownames(kingMat) <- king$sample.id

#look at the top corner of the matrix
kingMat[1:5,1:5]
kinship <- snpgdsIBDSelection(king)
head(kinship)

#run PC-AiR which provides robust population structure inference in samples with kinship and pedigree structure. 
pca <- pcair(gdsall, 
             kinobj = kingMat,
             kin.thresh=2^(-9/2),
             divobj = kingMat,
             div.thresh=-2^(-9/2))

#look at PCA results
names(pca)
length(pca$unrels)
head(pca$unrels)
head(pca$rels)

# extract the top 10 PCs and make a data.frame
pcs <- data.frame(pca$vectors[,1:10])
colnames(pcs) <- paste0('PC', 1:10)
pcs$sample.id <- pca$sample.id
dim(pcs)
head(pcs)

# PC-Relate, which provides accurate kinship inference even in the presence of population structure and ancestry admixture
seqData <- SeqVarData(gdsall)

# filter the GDS object to our LD-pruned variants
seqSetFilter(seqData, variant.id=pruned)

iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)

pcrel <- pcrelate(iterator, 
                  pcs=pca$vectors[,1:4], 
                  training.set=pca$unrels)
names(pcrel)

# relatedness between pairs of individuals
dim(pcrel$kinBtwn)
head(pcrel$kinBtwn)
# self-kinship estimates
dim(pcrel$kinSelf)
head(pcrel$kinSelf)

# transform the output into a kinship matrix .
pcrelMat <- pcrelateToMatrix(pcrel, scaleKin=1, verbose=FALSE)
dim(pcrelMat)


# load pheno data
phen <- read.table("Pheno_GWAS.txt", header = TRUE, sep = "\t", as.is = TRUE)
head(phen)

# select ids in this analysis
idsToKeep <- phen$sample.id
sampleIds <- data.frame(sample.id=seqGetData(gds, 'sample.id') )
phen <- dplyr::left_join(sampleIds, phen)

#Create an AnnotatedDataFrame including data and metadata. 
metadata <- data.frame(labelDescription = c(
  "subject identifier",
  "subject T2D",
  "subject's age in years",
  "subject's sex",
  "subject BMI"))

annot <- AnnotatedDataFrame(phen, metadata)
head(pData(annot))
varMetadata(annot)

class(annot$SEX)
annot$sex <- as.factor(annot$sex)
annot$T2D <- as.factor(annot$T2D)
annot$SEX <- as.numeric(as.character(annot$sex))
annot$T2D2 <- as.numeric(as.character(annot$T2D))
head(pData(annot))

# merge PCs with the sample annotation
dat <- left_join(pData(annot), pcs, by="sample.id")
# update the AnnotatedDataFrame
pData(annot) <- dat
head(pData(annot))

#fit regression model
nullmod <- fitNullModel(annot, 
                        outcome="T2D", 
                        covars=c("SEX", "age", "PC1", "PC2","PC3", "bmi"), 
                        cov.mat=kinship, 
                        family= "binomial")
nullmod$model
nullmod$fixef
head(nullmod$fit)

# run a single-variant test, accounting for genetic ancestry and genetic relatedness among the subjects. 
seqData <- SeqVarData(gds, sampleData=annot)
seqSetFilter(seqData, sample.id=idsToKeep)
iterator <- SeqVarBlockIterator(seqData, verbose=FALSE)
assoc <- assocTestSingle(iterator, nullmod,BPPARAM=BiocParallel::MulticoreParam())
dim(assoc)
head(assoc)

#add the alleles to the output
eff <- effectAllele(seqData, variant.id=assoc$variant.id)
assoc <- dplyr::left_join(assoc,eff)
head(assoc)

write.table(assoc,"/Output/GWAS_results.txt",row.names=FALSE,col.names=T,quote=FALSE,sep=" ")
nullmod$AIC
nullmod$logLik

#close the connection to the GDS file
seqClose(gdsall)

