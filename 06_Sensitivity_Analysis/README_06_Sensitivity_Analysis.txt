# Sensitivity Analysis
# Project: Epigenome-Wide Association Study of Adiponectin in West Africans

# Cohorts: RODAM-Pros (Ghanaians) and AADM (Nigerians)
# Author: Muhulo Muhau Mungamba
# Journal: eBioMedicine

#Overview

This folder contains all sensitivity analyses conducted to evaluate the robustness, interpretability, and validity of the main Epigenome-Wide Association Study (EWAS) findings reported in the manuscript.

All analyses in this folder are explicitly described in:

Methods → Sensitivity Analysis

Reviewer responses (Reviewers #1 and #2)

Each subfolder corresponds to a specific sensitivity question raised during peer review.

Folder Structure
06_Sensitivity_Analysis/
├── 06a_Geographical_Site_Specific_Variation/
├── 06b_Between_Cohort_Heterogeneity/
├── 06c1_Model_Collinearity_VIF_RODAM.R
├── 06c2_Model_Collinearity_VIF_AADM.R
└── README_06_Sensitivity_Analysis.txt

06a_Geographical_Site_Specific_Variation
Scientific question

Do adiponectin-associated CpG sites show differential methylation patterns across geographical environments?

Description

This analysis evaluates whether DNA methylation levels at genome-wide significant CpG sites vary across the three geographical locations represented in the RODAM cohort:

Rural Ghana

Urban Ghana

Amsterdam (Ghanaians)

Because AADM participants are all Nigerians without internal site stratification, this analysis is restricted to RODAM only, as stated in the manuscript.

Methods

DNA methylation beta values were extracted for top CpG sites.

Distributions were visualized using violin plots.

Kruskal–Wallis tests were used to assess global differences across sites.

Pairwise Wilcoxon tests (FDR-adjusted) were performed as a secondary step.

Scripts

06a1_ViolinPlots_KruskalWallis_RODAM.R

Output

Violin plots by site (Figure 6 in manuscript)

Kruskal–Wallis p-values

Pairwise comparison tables (Supplementary)

06b_Between_Cohort_Heterogeneity
Scientific question

Are the meta-analysis results consistent across cohorts?

Description

This analysis evaluates between-cohort heterogeneity for top CpG sites identified in the West African meta-analysis.

Methods

Heterogeneity statistics were extracted directly from METAL output.

Metrics include:

Cochran’s Q

I² statistic

Heterogeneity p-value

Scripts

06b_Heterogeneity_METAL.R

Output

Heterogeneity summary table for top CpG sites

Interpretation of heterogeneity levels (low / moderate / high)

Manuscript reference

Methods → Assessment of between-cohort heterogeneity

06c_Model_Collinearity
Scientific question

Does multicollinearity between BMI and adiponectin (or other covariates) bias EWAS results?

Description

This analysis addresses reviewer concerns regarding the inclusion of BMI as a covariate in the EWAS models, given its known biological relationship with adiponectin.

Analyses were performed separately for RODAM and AADM, as reported in the manuscript.

Components
1. BMI–Adiponectin correlation

Pearson correlation between BMI and log-transformed adiponectin

Scatter plots with linear regression lines

2️. Variance Inflation Factors (VIF)

VIFs computed for all covariates included in EWAS models

Adjusted GVIF values interpreted using standard thresholds:

<5: no collinearity

<10: acceptable

Scripts

06c1_Model_Collinearity_VIF_RODAM.R

06c2_Model_Collinearity_VIF_AADM.R

Manuscript reference

Methods → Model collinearity assessment

Reproducibility Notes

All analyses were performed in R (versions specified in individual scripts).

File paths in scripts reflect the original analysis environment and may require adaptation.

No individual-level genotype or methylation data are included in this repository.

Relationship to Main Analysis

These sensitivity analyses:

Do not redefine the main EWAS models

Do not alter the primary findings

Serve to validate modeling assumptions and support interpretability

They are included to ensure transparency and robustness in response to peer review.