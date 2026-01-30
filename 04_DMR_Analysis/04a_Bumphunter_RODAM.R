############################################
# Title: Differentially Methylated Regions (DMRs) using bumphunter
# Cohort: RODAM (Ghanaians)
# Outcome: Adiponectin EWAS
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-03-25
#
# Description:
# Identification of differentially methylated regions (DMRs)
# associated with log-transformed adiponectin levels using the
# bumphunter algorithm (minfi).
#
# Models adjusted for:
# Age, sex, BMI, recruitment site, batch, position, and estimated
# immune cell proportions.
#
# This analysis corresponds to the DMR section of the manuscript.
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

# Load required libraries
library(bumphunter)
library(GenomicRanges)

### ------------------------------------------------------------------
### Load normalized methylation data (RODAM)
### ------------------------------------------------------------------

# Normalized beta and M-values
beta <- read.table(
  "/project/Karlijn/EWAS_Adiponectin/Normalisation/20241025_Beta_Funnurm_NoSex_SNPs_NoCrosshybrRODAM_Adiponectin.txt",
  sep = "\t", header = TRUE
)

mval2 <- read.table(
  "/project/Karlijn/EWAS_Adiponectin/Normalisation/20241025_Mval_Funnurm_NoSex_SNPs_NoCrosshybrRODAM_Adiponectin.txt",
  sep = "\t", header = TRUE
)

# Probe annotation
probe.features <- read.table(
  "/project/Karlijn/EWAS_Adiponectin/Normalisation/20241025_Annot_NoXY_NoSNP_NoCrosshyb_Adiponectin_RODAM.txt",
  sep = "\t", header = TRUE
)

annotatie.gr <- makeGRangesFromDataFrame(
  probe.features,
  keep.extra.columns = TRUE,
  start.field = "pos",
  end.field = "pos"
)

### ------------------------------------------------------------------
### Load phenotype data
### ------------------------------------------------------------------

phenofile <- "/project/Karlijn/EWAS_Adiponectin/Pheno/20230105_RODAM_Adiponectin_Phenofile_wCells.csv"
targets <- read.table(phenofile, sep = ",", header = TRUE)

# Predictor
log_adip <- targets$log_Adiponectin

# Phenotypic covariates
sex  <- factor(targets$Sex)
age  <- targets$Age
BMI  <- targets$BMI
site <- factor(targets$Site)
DM   <- factor(targets$DM_Dichot)

# Technical covariates
batchArray <- factor(targets$BATCH_arry)
pos        <- targets$Sentrix_Position

# Estimated cell proportions
CD8T <- targets$CD8T
CD4T <- targets$CD4T
NK   <- targets$NK
Bcell<- targets$Bcell
Mono <- targets$Mono
Gran <- targets$Gran

### ------------------------------------------------------------------
### Design matrix (fully adjusted model)
### ------------------------------------------------------------------

design <- model.matrix(
  ~ log_adip + sex + age +
    CD8T + CD4T + NK + Bcell + Mono + Gran +
    batchArray + pos + BMI + site
)

### ------------------------------------------------------------------
### Define probe clusters
### ------------------------------------------------------------------

clust <- clusterMaker(
  pos = probe.features$pos,
  chr = probe.features$chr
)

### ------------------------------------------------------------------
### Step 1: Identify candidate regions (no permutation)
### ------------------------------------------------------------------

tab <- bumphunter(
  as.matrix(mval2),
  pos     = probe.features$pos,
  chr     = probe.features$chr,
  cluster = clust,
  design  = design,
  coef    = 2,
  cutoff  = 0.20,
  type    = "mval"
)

bumps1 <- tab$table
bumps1 <- subset(bumps1, L > 2)

### ------------------------------------------------------------------
### Step 2: Statistical inference with bootstrap
### ------------------------------------------------------------------

tab <- bumphunter(
  as.matrix(mval2),
  pos     = probe.features$pos,
  chr     = probe.features$chr,
  cluster = clust,
  design  = design,
  coef    = 2,
  cutoff  = 0.20,
  B       = 500,
  nullMethod = "bootstrap",
  type    = "mval"
)

bumps1 <- tab$table
bumps1 <- subset(bumps1, L > 2)

### ------------------------------------------------------------------
### Annotate DMRs with nearest genes
### ------------------------------------------------------------------

bumps.gr <- makeGRangesFromDataFrame(
  bumps1,
  keep.extra.columns = TRUE
)

nearest.genes <- annotatie.gr[nearest(bumps.gr, annotatie.gr)]

bumps.annotated <- cbind(
  bumps1,
  Gene_start = start(nearest.genes),
  Gene_end   = end(nearest.genes),
  Gene_Name  = nearest.genes$UCSC_RefGene_Name
)

### ------------------------------------------------------------------
### Save results
### ------------------------------------------------------------------

write.table(
  bumps.annotated,
  "2025.03.25_DMRs_Bumphunter_RODAM_Adiponectin_AgeSexCellTechBMIsite_n315.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

############################################
# END OF SCRIPT
############################################
