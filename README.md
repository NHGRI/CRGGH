# Epigenome-Wide Association Study of Adiponectin in West Africans
# Project Overview

This repository contains the full analysis pipeline supporting the manuscript:

"Epigenome-wide analysis identifies novel DNA methylation markers for circulating adiponectin" submitted to eBioMedicine

Authors
Muhulo Muhau Mungamba,2*a, Johanna Wijburg1a, Eva L. van der Linden1, Felix. P Chilunga1, Ayo P. Doumatey4, Amy R. Bentley4, Charles F. Hayfron-Benjamin1,3,5, Constance R. Sewani-Rusike2, Benedicta N. Nkeh-Chungag6, Rexford S. Ahima4, Charles Agyemang1,7, Peter Henneman8, Adebowale A. Adeyemo4, Charles N. Rotimi4, Karlijn A.C. Meeks1,4,9b

Affiliations
1.	Department of Public and Occupational Health, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam Public Health Research Institute, Amsterdam, The Netherlands.
2.	Department of Human Biology, Faculty of Medicine and Health Sciences, Walter Sisulu University, Mthatha, South Africa.
3.	Department of Vascular Medicine, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam Cardiovascular Sciences, Amsterdam, The Netherlands.
4.	Center for Research on Genomics and Global Health, National Human Genome Research Institute, National Institutes of Health, Bethesda, MD, USA.
5.	Departments of Physiology and Anaesthesia/Critical Care, University of Ghana Medical School, Korle Bu Teaching Hospital, Ghana.
6.	Department of Biological and Environmental Sciences, Faculty of Natural Sciences, Walter Sisulu University, Mthatha, South Africa.
7.	Division of Endocrinology, Diabetes and Metabolism, Department of Medicine, The Johns Hopkins University School of Medicine, Baltimore, MD, United States.
8.	Department of Human Genetics, Epigenetics, Amsterdam Reproduction and Development, research Institute, Amsterdam University Medical Centers, Amsterdam, The Netherlands.
9.	Division of Endocrinology, Diabetes and Nutrition, Department of Medicine, University of Maryland School of Medicine, Baltimore, MD, USA. <br>

<sup>*</sup>Joint authorship

<sup>a</sup>First corresponding author: <br>
* Muhulo Muhau Mungamba, Department of Public and Occupational Health, Amsterdam Public Health Research institute, Amsterdam University Medical Centers, University of Amsterdam, Amsterdam, The Netherlands. Email: m.m.mungamba@amsterdamumc.nl / mmungamba@wsu.ac.za.

<sup>b</sup>Second corresponding author <br>
* Karlijn A.C. Meeks, Division of Endocrinology, Diabetes and Nutrition, Department of Medicine, University of Maryland School of Medicine, 670 W. Baltimore Street, HSF III 4061, Baltimore, MD 21201, USA. Email: karlijn.meeks@som.umaryland.edu 
 
# Abstract

Background: 
Adiponectin is a circulating adipokine involved in energy metabolism and inflammation, with reported protective effects against cardiometabolic diseases such as type 2 diabetes (T2D) and early kidney disease. However, its regulation remains poorly understood. This study aimed to identify epigenetic loci associated with adiponectin levels.

Methods: 
DNA methylation was profiled using Illumina 450K and EPIC (850K) arrays in 315 Ghanaians (RODAM-Pros study) and 593 Nigerians (AADM study). Differentially methylated positions (DMPs) were identified using linear regression models adjusted for age, sex, BMI, blood cell proportions, and technical covariates. Analyses were stratified by T2D status and cohort, then meta-analyzed to identify DMPs associated with adiponectin across T2D status (combining participants with and without diabetes). RNA-seq data on 77 blood, 49 subcutaneous adipose tissue (SAT), and 55 skeletal muscle samples from the AADM study were used to identify eQTMs for identified DMPs.

Findings: 
We identified three epigenome-wide significant DMPs: cg03546163 (Z-score=5.76, p=<0.001, 5’ UTR of FKBP5), cg02561343 (Z-score=5.11, p=<0.001, within UST), and cg23969380 (Z-score=5.13, p<=0.001, ADGRD1 body). cg03546163 was an eQTM for PLA2G12B in SAT (beta=-0.039, FDR=0.047), cg02561343 for PSMD8 (beta=-11.85, FDR=0.029) and TECR (beta=-9.48, FDR=0.029) in SAT, and cg23969380 for HIGD2AP1 (beta=-0.095, FDR=0.024) in blood. These genes have been reported to be involved in lipid metabolism (PLA2G12B and TECR), proteasomal degradation (PSMD8), and cellular-stress responses (HIGD2AP1).

Interpretation: 
This first epigenome-wide study of adiponectin in sub-Saharan African populations identified novel DNA methylation loci potentially involved in adiponectin regulation through lipid metabolism, inflammation, proteostasis, and stress response pathways. These findings provide a foundation for replication and further investigation to understand the role of adiponectin in cardiometabolic health.  

# Study Cohorts

RODAM-Pros: Ghanaian migrants and non-migrants (rural Ghana, urban Ghana, Amsterdam)

AADM: Nigerian adults

Analyses were stratified by type 2 diabetes (T2D) status where appropriate, reflecting the strong biological relationship between adiponectin, adiposity, and T2D.

# Software Environment

All analyses were performed in the R statistical computing environment.

# Key R packages include:

minfi, limma, bacon, missMethyl, clusterProfiler, dmrff, bumphunter,
enrichR, car, ggplot2, dplyr

Meta-analyses were performed using METAL (Linux command line).

# Repository Structure
├── 01_RODAM/                 # RODAM cohort: phenotype prep, normalization, EWAS <br>
├── 02_AADM/                  # AADM cohort: EWAS documentation & inputs <br>
├── 03_MetaAnalysis/          # METAL meta-analyses + QC plots <br>
├── 04_DMR_Analysis/          # Region-level analyses (bumphunter, DMRff) <br>
├── 05_InSilico_Analysis/     # Functional annotation & enrichment <br>
├── 06_Sensitivity_Analysis/  # Robustness & reviewer-requested analyses <br>
└── README.txt                # This file <br>


Each folder contains a dedicated README describing:
* the scientific question addressed
* methods used
* outputs generated
* how the analyses map directly to the manuscript

#### Folder Descriptions

# 01_RODAM

Complete analysis pipeline for the RODAM-Pros cohort, including:
* Phenotype preparation
* DNA methylation quality control
* Functional normalization
* Probe filtering (sex chromosomes, SNPs, cross-hybridization)
* EWAS stratified by T2D status (cases, controls, combined)
* Fully adjusted models including age, sex, BMI, cell counts, and technical covariates
* QC visualizations (QQ plots, Manhattan plots, genomic inflation)
This folder is fully self-contained and reproducible with access to RODAM data.

# 02_AADM

Documentation and inputs for the AADM cohort EWAS. <br>
EWAS was conducted by a collaborating analyst (Karlijn A.C. Meeks)<br>
Raw EWAS scripts are not duplicated here<br>
Summary statistics and metadata required for downstream analyses are included<br>
All subsequent meta-analysis, in-silico, and sensitivity analyses rely on these outputs<br>
Note: The absence of raw AADM EWAS scripts reflects a collaborative workflow and does not affect reproducibility of the reported results.<br>

# 03_MetaAnalysis

Fixed-effects inverse-variance weighted meta-analyses performed using METAL.
Analyses include:
* West Africans (RODAM + AADM; cases + controls) — main results
* T2D cases only
* T2D controls only
* Cohort-specific summaries (RODAM only, AADM only)
* Outputs include:
* QQ plots and genomic inflation (lambda)
* Manhattan plots (BACON-corrected p-values)
* Between-cohort heterogeneity statistics (Q, I²)

# 04_DMR_Analysis

Region-level methylation analyses using:
* bumphunter (cohort-specific DMRs)
* DMRff (meta-analysis summary statistics)
Includes:
* Genomic annotation of DMRs
* Region characterization
* Methods aligned with the manuscript and reviewer requests

# 05_InSilico_Analysis

Functional interpretation of EWAS findings, including:
* Pathway enrichment (missMethyl, clusterProfiler)
* Transcription factor enrichment (enrichR)
* Public eQTM resources (iMETHYL)
* Regulatory feature enrichment using eFORGE (web-based; documented)
* eQTM analyses using individual-level expression data were performed externally and are documented accordingly.

# 06_Sensitivity_Analysis

Sensitivity analyses conducted to evaluate robustness and address reviewer concerns, including:
* Geographical site-specific variation (RODAM only)
* Between-cohort heterogeneity
* Model collinearity assessment (BMI and covariates)
All analyses are explicitly referenced in the manuscript and reviewer response letter. <br>
See README_06_Sensitivity_Analysis.txt for detailed descriptions.

# Reproducibility Notes

* Scripts reflect the original computational environments and directory structures <br>
* File paths may require modification before reuse <br>
* Individual-level methylation and phenotype data are not publicly shared due to ethical and consent restrictions
* All results reported in the manuscript can be regenerated using the provided scripts with appropriate data access

# Relationship to the Manuscript

This repository corresponds directly to the following manuscript sections:
* Statistical Analysis
* Differentially Methylated Positions (DMPs)
* Differentially Methylated Regions (DMRs)
* In-silico Functional Annotation
* Sensitivity Analyses
Scripts are organized to mirror the manuscript structure and reviewer-requested additions.

# Contact (Corresponding Authors)

Corresponding author:<br>
Muhulo Muhau Mungamba<br>
Department of Public and Occupational Health<br>
Amsterdam Public Health Research Institute<br>
Amsterdam University Medical Centers, University of Amsterdam<br>
Email: m.m.mungamba@amsterdamumc.nl/mmungamba@wsu.ac.za

Second corresponding author<br>
Karlijn A.C. Meeks<br>
Division of Endocrinology, Diabetes and Nutrition<br>
Department of Medicine<br><br>
University of Maryland School of Medicine
Baltimore, MD, USA<br>
Email: karlijn.meeks@som.umaryland.edu<br>




