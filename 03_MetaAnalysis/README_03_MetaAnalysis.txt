# Meta-analysis of Adiponectin EWAS

## Project
Meta-analysis of EWAS results from RODAM-Pros and AADM

## Author
Muhulo Muhau Mungamba

## Overview

This folder contains all scripts used to perform fixed-effects meta-analyses of adiponectin-associated DNA methylation across cohorts.

Meta-analyses were conducted using METAL (Linux command line) and followed by downstream QC and visualization in R.

## Meta-analyses Performed

### 1. West Africans (Main analysis)
- RODAM + AADM
- T2D cases + controls combined
- Represents the primary results reported in the manuscript

### 2. Stratified meta-analyses
- T2D cases only (RODAM + AADM)
- T2D controls only (RODAM + AADM)

### 3. Cohort-specific summaries
- RODAM only (cases + controls)
- AADM only (cases + controls)

## Statistical Methods

- Fixed-effects inverse-variance weighted meta-analysis
- Implemented using METAL
- Justification:
  - Two cohorts only
  - Random-effects models would be unstable
- Between-cohort heterogeneity assessed using:
  - Cochran’s Q
  - I² statistic

## Post-meta-analysis QC

- QQ plots and genomic inflation factors (λ)
- Manhattan plots (BACON-corrected p-values)
- Separate plots generated for each meta-analysis

## Manuscript Reference

Methods → Meta-analysis  
Results → West African meta-analysis findings  
Supplementary Figures (QQ plots)

## Notes

- Meta-analysis performed on 409,033 shared CpG sites
- Genomic build: GRCh37 / hg19
- METAL scripts were executed in Linux; R scripts handle visualization and interpretation

#############################################################################################################################################################################################################################################