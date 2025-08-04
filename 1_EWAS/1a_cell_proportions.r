##############################
##Calculate Cell Proportions##
##############################

## load packages
library(IlluminaHumanMethylationEPICmanifest)
library(minfi)
library("FlowSorted.Blood.EPIC")

### set the work directory
workdirectory="/Cells"
setwd(workdirectory)

options(stringsAsFactors = F)

### read phenotype file
targets = read.csv("/Pheno/Phenotype_file.csv", header=T)
head(targets)


### read the idat data of the samples in the phenotype file
idatspath="/idats" 
BasePath=idatspath
RGset=read.metharray.exp(BasePath,targets) ## read datafiles
pd <- pData(RGset) ## accesses phenoData

dim(RGset)

### Calculate the blood distributions
x1=estimateCellCounts2(RGset, compositeCellType = "Blood",
                      cellTypes = c("CD8T","CD4T", "NK","Bcell","Mono","Neu"), referencePlatform = "IlluminaHumanMethylationEPIC",
                      returnAll = T, meanPlot = F, verbose = F)

## Extract table of proportions
x2 <- print(x1$prop)


## Add "Basename" as first column header
x2_df <- as.data.frame(x2) #make it into a dataframe first
basenames <- rownames(x2_df) #rownames are the basenames
x2_df$Basename <- basenames
x2_df <- x2_df[, c("Basename", colnames(x2_df)[-ncol(x2_df)])] #Reorder the columns
head(x2_df)


## Write cell table to csv
write.table(x2_df,"cell_proportions.csv",sep=",", quote=F, row.names=FALSE)


## Merge cell proportions with pheno file for EWAS
pheno_celldistr = merge(targets,x2_df,by=c("Basename"),all.x=T)


## Write pheno with cells to CSV file
write.csv(pheno_celldistr, "Pheno_wCells.csv", quote = F, row.names=F) 







