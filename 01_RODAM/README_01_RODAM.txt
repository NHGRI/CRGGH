# RODAM EWAS – Adiponectin

## Project
Epigenome-Wide Association Study (EWAS) of circulating adiponectin levels  
Cohort: RODAM-Pros (Ghanaians)

## Author
Muhulo Muhau Mungamba

## Overview

This folder contains all scripts used to perform the cohort-specific EWAS of adiponectin in the RODAM-Pros study. These analyses form a core component of the manuscript and include phenotype preparation, normalization, EWAS modeling, and visualization.

All analyses are explicitly described in the Methods section of the manuscript.

## Analyses Included

### 1. Phenotype preparation
- Cleaning and transformation of adiponectin levels
- Log-transformation of adiponectin
- Integration of covariates:
  - Age, sex, BMI
  - Recruitment site
  - Technical variables (batch, plate position)
  - Estimated immune cell proportions

### 2. DNA methylation preprocessing
- Functional normalization (minfi)
- Probe filtering:
  - Sex chromosomes
  - SNP-associated probes
  - Cross-hybridizing probes
- Generation of normalized M-values and beta-values

### 3. EWAS modeling
- Linear models fitted using limma
- DNA methylation (M-values) as outcome
- Log-transformed adiponectin as exposure
- Fully adjusted models including biological and technical covariates
- Stratification by T2D status:
  - T2D cases
  - T2D controls
  - Combined (cases + controls)

### 4. Quality control and visualization
- QQ plots and genomic inflation (λ)
- Bacon correction for test statistic inflation
- Manhattan plots (BACON-corrected p-values)

## Output

- EWAS summary statistics (TopTables)
- QQ plots and lambda values
- Manhattan plots used in the manuscript and supplementary figures

## Manuscript Reference

Methods → Differentially methylated positions (DMPs)  
Results → Cohort-specific EWAS findings

## Notes

- M-values were used for statistical inference.
- Beta-values were used for visualization only.
- All file paths reflect the original analysis environment and may require adaptation.

#########################################################################################################################################################################