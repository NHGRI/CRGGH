# AADM EWAS – Adiponectin

## Project
Epigenome-Wide Association Study (EWAS) of circulating adiponectin levels  
Cohort: AADM (Nigerians)

## Overview

This folder documents the AADM cohort-specific EWAS analyses used in the manuscript.

The EWAS for AADM was performed by a collaborating analyst (Karlijn Meeks).  
Raw EWAS scripts are therefore not reproduced here.

Downstream analyses (meta-analysis, in-silico annotation, sensitivity analyses) rely on the resulting EWAS summary statistics, which are fully documented and reproducible.

## EWAS Description (AADM)

- Linear EWAS models fitted using limma
- DNA methylation M-values as outcome
- Log-transformed adiponectin as exposure
- Covariate adjustment:
  - Age
  - Sex
  - BMI
  - Estimated immune cell proportions
- Stratification by T2D status:
  - T2D cases
  - T2D controls
  - Combined

These models are described in detail in the manuscript Methods.

## Data Used in This Repository

The following AADM EWAS outputs are used in later folders:

- EWAS summary statistics (TopTables)
- Effect estimates (β), standard errors, and p-values
- Shared CpGs between EPIC (AADM) and 450K (RODAM)

## Manuscript Reference

Methods → Differentially methylated positions (DMPs)  
Meta-analysis section

## Reproducibility Note

The absence of raw AADM EWAS scripts reflects the collaborative nature of the project and does not affect reproducibility of the reported results, as all downstream analyses are based on shared summary statistics.

##########################################################################################################################################################################################################################