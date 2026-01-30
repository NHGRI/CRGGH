# Differentially Methylated Region (DMR) Analysis

## Project
Region-level DNA methylation analysis for adiponectin EWAS

## Author
Muhulo Muhau Mungamba

## Overview

This folder contains scripts used to identify differentially methylated regions (DMRs) associated with adiponectin levels.

Both cohort-specific and meta-analysis-based region detection approaches were applied, as described in the manuscript.

## Analyses Included

### 1. Cohort-specific DMRs (bumphunter)

- Performed separately in:
  - RODAM-Pros
  - AADM
- Implemented using the bumphunter function (minfi)
- Regions defined as:
  - ≥3 adjacent CpGs
  - Effect size cutoff applied
- Statistical significance:
  - Family-wise error rate (FWER) < 0.25

### 2. Meta-analysis DMRs (DMRff)

- Applied to METAL summary statistics
- Enables region detection without individual-level data
- Adjacent CpGs (<500 bp) grouped into regions
- FDR correction applied
- Regions annotated using UCSC (GRCh37/hg19)

## Outputs

- Cohort-specific DMR tables
- Meta-analysis DMR tables
- Annotated DMR summary tables
- Supplementary tables reported in the manuscript

## Manuscript Reference

Methods → Differentially methylated regions (DMRs)  
Results → DMR findings

## Notes

- bumphunter analyses use individual-level methylation data
- DMRff analyses use meta-analysis summary statistics
- Visualization of DMRs was limited to descriptive plots and summary tables

###########################################################################################################################################################################################################################