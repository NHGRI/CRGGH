###############
##MR Analyses##
###############
  
# load the libraries
library(here)
library(MendelianRandomization)
library(knitr)
library(readr)
library(tidyverse)
library(tableone)
library(dplyr)
library(ggplot2)

# load extracted instrument snps
extr <- read_table(here("Output/Extracted_variants", "mQTL_FDR0.05.r01.clumped.extracted.txt"))
nrow(extr)

# read pheno data
df1 <- read_csv(here("Pheno", "Pheno.csv"))
df1 <- as_tibble(df1)
df1$sample.id <- df1$vcf_id
nrow(df1)
dim(df1)
head(df1)

# select pheno variables to test
myvars <- c("sample.id","age","sex", "waist", "hip", "bmi", "smoking")
df1x <- df1[myvars]

# add pheno to variants
dfextr <- merge(df1x,extr,by=c("sample.id"))
head(dfextr)
nrow(dfextr)

# factors
dfextr$sex <-as.factor(dfextr$sex)
dfextr$smoking <-as.factor(dfextr$smoking)

## Associate instruments with confounders

# Initialize an empty dataframe to store p-values
p_values_df <- data.frame()

# Loop over the response variables (columns 9 to the last column)
for (i in 9:ncol(dfextr)) {
  # Temporary list to store p-values for the current response variable
  temp_p_values <- c()
  
  # Loop over the predictor variables (columns 3 to 8)
  for (j in 2:7) {
    # Create the formula dynamically with one predictor at a time
    formula <- as.formula(paste(colnames(dfextr)[i], "~", colnames(dfextr)[j]))
    
    # Perform the linear regression
    model <- lm(formula, data = dfextr)
    
    # Extract the p-value for the predictor
    p_value <- summary(model)$coefficients[2, 4]  # 2nd row corresponds to the predictor variable
    temp_p_values <- c(temp_p_values, p_value)
  }
  
  # Add the p-values for the current response variable to the dataframe
  p_values_df <- rbind(p_values_df, temp_p_values)
}

# Set the column names to the predictor variables
colnames(p_values_df) <- colnames(dfextr)[2:7]

# Set the row names to the response variables
rownames(p_values_df) <- colnames(dfextr)[9:ncol(dfextr)]

# View the p-values dataframe with row and column names
head(p_values_df)


## Extract violating instruments

# make a list of instruments associated with any of the confounders at p<0.05
sig5 <- subset(p_values_df, p_values_df$age < 0.05 | p_values_df$sex < 0.05 | p_values_df$bmi < 0.05 | p_values_df$waist < 0.05 | p_values_df$hip < 0.05 | p_values_df$smoking < 0.05)
nrow(sig5)
sig5 <- rownames_to_column(sig5, var = "variant")
print(sig5)


## Getting MR input data

# load CIS meQTLs
mqtl <- read_table(here("Output/mQTL", "mQTLs_out_fdr.05_57CpGs.txt"))
nrow(mqtl)
head(mqtl)

# load the clumps from clumping for T2D
clump <- read_table(here("/Clumping", "allcpgs_mqtls_fdr0.05.r01.clumped.txt"))
nrow(clump)
clump$variant <- clump$SNP
head(clump)
nrow(clump)

# load SNP - outcome associations
gen <- read_table(here("Output/GWAS", "GWAS_results.txt" ))
head(gen)
nrow(gen)
gen$variant <- gen$variant.id


## Merge files

# add mQTL data to only the clumps
myvars <- c("variant")
clump2 <- clump[myvars]

dfmcl <- merge(clump2,mqtl,by=c("variant"), all.x=TRUE)
head(dfmcl)
nrow(dfmcl)

dfgen <- merge(dfmcl,gen,by=c("variant"))
head(dfgen)
nrow(dfgen)


## Exclude violating instruments

sig5 <- sig5 %>%
  mutate_all(~ gsub("_", ":", .))
del <- as.list(sig5$variant)
del

dfgen2 <- dfgen %>%
  filter(!variant %in% del)
nrow(dfgen2)
dfgen <- dfgen2

# Make sure variables are numeric
dfgen$beta <- as.numeric(dfgen$beta)
dfgen$beta_se <- as.numeric(dfgen$beta_se)

## Minor allele check

#QTL tools uses the ALT allele from the VCF as the effect allele.

#Make sure that the effect allele from SNP-outcome associations is the same as the minor allele from the SNP-exposure associations

# Separate the 'variant' column into 'chr', 'pos', 'alt', and 'ref'
dfgen <- dfgen %>%
  separate(variant, into = c("chr", "pos", "alt", "ref"), sep = ":", remove = FALSE)

#If effect allele is not the same as minor, multiply GWAS estimate by -1
dfgen <- dfgen %>%
  mutate(Est = if_else(alt != effect.allele, Est * -1, Est))


## Align direction of effects

#make sure that all snp-exposure associations have the same direction of effect (variant associated with higher methylation)
ref_direction <- ifelse(dfgen$beta > 0, 1, -1)
dfgen$beta <- dfgen$beta * ref_direction
dfgen$se <- dfgen$se * abs(ref_direction)

#align the snp-outcome associations
dfgen$Est <- dfgen$Est * ref_direction
dfgen$Est.SE <- dfgen$Est.SE * abs(ref_direction)



## Create tables for MR analyses
cpgs <- unique(dfgen$cpg)

unique_cpg_values <- unique(dfgen$CpG)

for (cpg_value in unique_cpg_values) {
  df_filtered <- filter(dfgen, CpG == cpg_value)
  assign(paste("df", cpg_value, sep = "_"), df_filtered)
}

## Fixed Effect Model

#Step 1 of the Rucker framework is the fixed effect  model

results_list <- list()

# Loop through each unique CpG value
for (cpg_value in unique_cpg_values) {
  
  # Generate the dataframe name
  df_name <- paste("df", cpg_value, sep = "_")
  
  # Create mrobj object
  mrobj_name <- paste("mrobj", cpg_value, sep = "_")
  
  # Access the dataframe using the generated name
  df <- get(df_name)
  
  # Create mrobj using the dataframe columns
  assign(mrobj_name, mr_input(bx = df$beta, bxse = df$se,
                              by = df$Est, byse = df$Est.SE))
  
  # Perform mr_ivw using the created mrobj
  result <- mr_ivw(get(mrobj_name), model = "fixed")
  
  # Extract required values from the result
  n_variants <- result$SNPs
  fstat <- result$Fstat
  estimate <- result$Estimate
  std_error <- result$StdError
  ci_lower <- result$CILower
  ci_upper <- result$CIUpper
  p_value <- result$Pvalue
  Q <- result$Heter.Stat[1]  # Cochran's Q statistic
  Q_pval <- result$Heter.Stat[2]  # Cochran's Q p-value
  
  # Store the values in the list as a named vector
  results_list[[cpg_value]] <- c(
    n_variants = n_variants,
    F_statistic = fstat,
    Estimate = estimate,
    StdError = std_error,
    CILower = ci_lower,
    CIUpper = ci_upper,
    Pvalue = p_value,
    Cochran_Q = Q,
    Cochran_Q_Pval = Q_pval
  )
}

# Convert the list to a dataframe
results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x))))
results_df <- data.frame(CpG = names(results_list), results_df)

# Check the dataframe
head(results_df)



## Random Effect Model

#Random effects model is step 2 in the Rucker framework
#only needed for CpGs with substantial heterogeneity in step 1

results_list <- list()

# Loop through each unique CpG value
for (cpg_value in unique_cpg_values) {
  
  # Generate the dataframe name
  df_name <- paste("df", cpg_value, sep = "_")
  
  # Create mrobj object
  mrobj_name <- paste("mrobj", cpg_value, sep = "_")
  
  # Access the dataframe using the generated name
  df <- get(df_name)
  
  # Create mrobj using the dataframe columns
  assign(mrobj_name, mr_input(bx = df$beta, bxse = df$se,
                              by = df$Est, byse = df$Est.SE))
  
  # Perform mr_ivw using the created mrobj
  result <- mr_ivw(get(mrobj_name), model = "random")
  
  # Extract required values from the result
  n_variants <- result$SNPs
  fstat <- result$Fstat
  estimate <- result$Estimate
  std_error <- result$StdError
  ci_lower <- result$CILower
  ci_upper <- result$CIUpper
  p_value <- result$Pvalue
  Q <- result$Heter.Stat[1]  # Cochran's Q statistic
  Q_pval <- result$Heter.Stat[2]  # Cochran's Q p-value
  
  # Store the values in the list as a named vector
  results_list[[cpg_value]] <- c(
    n_variants = n_variants,
    F_statistic = fstat,
    Estimate = estimate,
    StdError = std_error,
    CILower = ci_lower,
    CIUpper = ci_upper,
    Pvalue = p_value,
    Cochran_Q = Q,
    Cochran_Q_Pval = Q_pval
  )
}

# Convert the list to a dataframe
results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x))))
results_df <- data.frame(CpG = names(results_list), results_df)

# Check the dataframe
head(results_df)




## MR Egger

#MR Egger only for CpGs with large Rucker's Q statistic

results_list <- list()

# Loop through each unique CpG value
for (cpg_value in unique_cpg_values) {
  
  # Generate the dataframe name
  df_name <- paste("df", cpg_value, sep = "_")
  
  # Create mrobj object
  mrobj_name <- paste("mrobj", cpg_value, sep = "_")
  
  # Access the dataframe using the generated name
  df <- get(df_name)
  
  # Create mrobj using the dataframe columns
  assign(mrobj_name, mr_input(bx = df$beta, bxse = df$se,
                              by = df$Est, byse = df$Est.SE))
  
  # Perform mr_ivw using the created mrobj
  result <- mr_egger(get(mrobj_name))
  
  # Extract required values from the result
  n_variants <- result$SNPs
  estimate <- result$Estimate
  std_error <- result$StdError.Est
  ci_lower <- result$CILower.Est
  ci_upper <- result$CIUpper.Est
  p_value <- result$Pvalue.Est
  intercept <- result$Intercept  # Mr Egger intercept
  interc_pval <- result$Pvalue.Int  # Intercept p-value
  Q <- result$Heter.Stat[1]  # Cochran's Q statistic
  Q_pval <- result$Heter.Stat[2]  # Cochran's Q p-value
  
  # Store the values in the list as a named vector
  results_list[[cpg_value]] <- c(
    n_variants = n_variants,
    Estimate = estimate,
    StdError = std_error,
    CILower = ci_lower,
    CIUpper = ci_upper,
    Pvalue = p_value,
    Intercept = intercept,
    Intercept_pval = interc_pval,
    Cochran_Q = Q,
    Cochran_Q_Pval = Q_pval
  )
}

# Convert the list to a dataframe
results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x))))
results_df <- data.frame(CpG = names(results_list), results_df)

# Check the dataframe
head(results_df)



## Rucker Q statistic

# run ivw random and mr egger for the two cpgs with significant Q-statistic

###cg26872670
ivw_cg26872670 <- mr_ivw(mrobj_cg26872670, model = "random")
egger_cg26872670 <- mr_egger(mrobj_cg26872670)

# Extract Q IVW and Q Egger
Q_ivw <- ivw_cg26872670$Heter.Stat[1]
Q_egger <- egger_cg26872670$Heter.Stat[1]

# Calculate Rucker's Q statistic
Q_Rucker <- Q_ivw - Q_egger

# Print the result
cat("Rucker's Q statistic cg26872670:", Q_Rucker)


###cg26591384
ivw_cg26591384 <- mr_ivw(mrobj_cg26591384, model = "random")
egger_cg26591384 <- mr_egger(mrobj_cg26591384)

# Extract Q IVW and Q Egger
Q_ivw <- ivw_cg26591384$Heter.Stat[1]
Q_egger <- egger_cg26591384$Heter.Stat[1]

# Calculate Rucker's Q statistic
Q_Rucker <- Q_ivw - Q_egger

# Print the result
cat("Rucker's Q statistic cg26591384:", Q_Rucker)

