#################
##Normalisation##
#################

## load packages
library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(minfiDataEPIC)
library("FlowSorted.Blood.EPIC")
library(affy)
library(limma)

options(stringsAsFactors = F)


### set working directory
workdirectory="/Data" 
setwd(workdirectory)

## load idat path
idatspath="/idats" 
BasePath=idatspath

## load pheno with calculated cell proportions
pheno="/Pheno/Pheno_wCells.csv"

## load cross-hybridization probes
cross="/Data/CrosHYB_SNPSlist.csv"

## read phenotypefile
targets <- read.table(pheno, sep=",",header=T)
head(targets)

## read the idat data of the samples in the phenotype file
RGset=read.metharray.exp(BasePath,targets,force=T)

### get unnormalized beta values
betaraw=getBeta(RGset) 

#### normalize data using Funnorm
GMset=preprocessFunnorm(RGset, nPCs=2, sex = c("male","female"), bgCorr = TRUE,
                        dyeCorr = TRUE, verbose = TRUE)

pd <- pData(RGset) 


## get annotation to filter probes
annotation <- getAnnotation(GMset)

## make a matrix from the annotation
probe.features=as.matrix(annotation)

## Get X and Y Probes
indices <- which(annotation$chr == "chrX"| annotation$chr == "chrY")

## filter out X and Y probes
GMset.nosex <- GMset[-indices,] #Probes remaining n= 846,232

## filterout probes with snps & cross hybridization
probes4removal <- read.csv(cross, stringsAsFactors  = F, col.names = 1)[,1]
GMset.culled <- GMset.nosex [which(!featureNames(GMset.nosex) %in% probes4removal), ] #Probes remaining n= 803,022


## Prepare Betas and M-values 
mval = getM(GMset.culled)
beta2=getBeta(GMset.culled)
x=which(mval == -Inf, arr.ind = T)
xx=sort(rownames(x))
g=mval
zz1=g[!rownames(g) %in% xx, ]
mval2=zz1
beta=beta2[!rownames(beta2) %in% xx, ]


## write betas and M values to file
write.table(mval2,"Mvalues.txt", sep="\t")
write.table(beta,"Betas.txt", sep="\t")

## write probe annotation to file, excl. XY, noSNP, NoCrosshyb
probe.features=as.matrix(annotation)
x=data.frame(probe.features[match(row.names(beta), row.names(probe.features)),])## get annotation data for probes
write.table(x,"Annotation.txt", sep="\t")


