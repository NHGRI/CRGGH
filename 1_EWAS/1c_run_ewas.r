############
##Run EWAS##
############

## load libraries

library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library(minfiDataEPIC)
library(limma)
library(affy)

options(stringsAsFactors = F)

### set working directory
workdirectory="/Output"  

## load idat path
idatspath="/idats" 

## Read normalised Beta and M-values from step 2
beta <- read.table("/Data/Betas.txt", sep="\t",header=T)
mval2 <- read.table("/Data/Mvalues.txt", sep="\t",header=T)
probe.features <- read.table("/Data/Annotation.txt", sep="\t",header=T)

## load pheno with calculated cell proportions
pheno="/Pheno/Pheno_wCells.csv"

targets <- read.table(phenofile, sep=",",header=T)
head(targets)
 
## prepare covariates
T2D = factor(targets$T2D)
HBA1C = targets$log_hba1c
HOMAS = targets$log_homas
HOMAB = targets$log_homab
sex=factor(targets$sex)
age=targets$age
site=factor(targets$site)
cd8t=targets$CD8T
cd4t=targets$CD4T
nk=targets$NK
bcell=targets$Bcell
mono=targets$Mono
neu=targets$Neu
BMI=targets$bmi
pos=targets$Sentrix_Position
plate=factor (targets$Sample_Plate)


## Design model

### binary traits (i.e. T2D)
design= model.matrix(~0+T2D+sex+age+BMI+cd8t+cd4t+nk+bcell+mono+neu+batch+pos+plate)
design1= model.matrix(~0+T2D) ### model matrix without covariates for comparison

cont.matrix=makeContrasts(caseContr=(group1-group0), levels=design)# make contrasts, 1 is case 0 is control
cont.matrix2=makeContrasts(caseContr=(group1-group0), levels=design1)

### continuous traits (i.e. HBA1c, HOMAS, and HOMAB)
design= model.matrix(~homas+sex+age+BMI+cd8t+cd4t+nk+bcell+mono+neu+plate+pos)
design1= model.matrix(~homas)


## analyses with M-values
fit=lmFit(mval2, design)
fitN=lmFit(mval2, design1)

fit2=contrasts.fit(fit,cont.matrix) #for binary models
fit2N=contrasts.fit(fitN,cont.matrix2) #for binary models

fit2e=eBayes(fit2)
fit2eN=eBayes(fit2N)

topTableFDR= topTable(fit2e, number=nrow(mval2), adjust="fdr", confint=TRUE)
topTable= topTable(fit2eN, number=nrow(mval2),adjust="none", confint=TRUE)# design1= model.matrix(~syst)


## Add annotation
x=data.frame(probe.features[match(row.names(topTableFDR), row.names(probe.features)),])## get annotation data for probes

z=data.frame(topTable[match(row.names(topTableFDR),row.names(topTable)),])

y=z$logFC
topTableFDRannot= cbind(topTableFDR,y)

colnames(topTableFDRannot)[9]="DeltaMval"
colnames(topTableFDRannot)[1]="logFCMval"

## Add SE
topTableFDRannot$SE_mval=(topTableFDRannot$CI.R-topTableFDRannot$CI.L)/3.92

## write toptable to a file
write.table(topTableFDRannot,"/Output/TopTable.txt", quote = FALSE, sep="\t",col.names=T,row.names=T) ## write data to file



