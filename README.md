<img align="right" width="200" height="200" src="https://github.com/user-attachments/assets/5361ef5a-1540-4e73-b555-75a4441ebb3c"> <br>

# Causal relationship between epigenetic markers and type 2 diabetes in West Africans: A Mendelian randomisation analysis

##### Authors:
###### Karlijn A.C. Meeks(1,2,3,4), Eva L. van der Linden(4), Amy R. Bentley(1), Ayo P. Doumatey(1), Peter Henneman(5), Nora Franceschini(6), Themistocles L. Assimes(7,8), Felix P. Chilunga(4), Charles F. Hayfron-Benjamin(4,9), Ellis Owusu-Dabo(10), Guanjie Chen(1), Charles Agyemang(4), Adebowale A. Adeyemo(1), Charles N. Rotimi(1)<br>
###### 1. Center for Research on Genomics and Global Health, National Human Genome Research Institute, National Institutes of Health, Bethesda, MD, 20892, USA
###### 2. Division of Endocrinology, Diabetes and Nutrition, Department of Medicine, University of Maryland School of Medicine, Baltimore, MD, USA
###### 3. Department of Epidemiology and Public Health, University of Maryland School of Medicine, Baltimore, MD, USA
###### 4. Department of Public and Occupational Health, Amsterdam UMC, University of Amsterdam, Amsterdam Public Health Research Institute, Amsterdam, The Netherlands
###### 5. Department of Human Genetics, Epigenetics, Reproduction & Development, Amsterdam UMC, University of Amsterdam, Amsterdam, The Netherlands
###### 6. Departments of Epidemiology and Genetics, University of North Carolina, Chapel Hill, NC, USA
###### 7.	Department of Medicine, Stanford University School of Medicine, Stanford, CA, USA
###### 8.	VA Palo Alto Healthcare System, Palo Alto, CA, USA
###### 9.	Department of Physiology, University of Ghana Medical School, Accra, Ghana
###### 10.	Department of Global and International Health, School of Public Health, Kwame Nkrumah University of Science and Technology, Kumasi, Ghana

## Abstract <br>
Aims: Evidence for a causal role of DNA methylation sites (CpGs) in type 2 diabetes (T2D) and glycaemic traits is limited due to the cross-sectional nature of many epigenome-wide association studies (EWAS). In addition, epigenetic studies in West African populations are particularly sparse despite the high and rising burden of T2D in these populations. Hence, we aimed to identify CpGs causally associated with T2D among West Africans by leveraging Mendelian randomisation (MR) analysis and longitudinal data. 

Methods: We used the Illumina EPIC DNA-methylation array to profile the methylation of DNA extracted from white blood cells collected from 879 Ghanaians (Research on Obesity and Diabetes among African Migrants study [RODAM]) and 332 Nigerians (Africa America Diabetes Mellitus study [AADM]) who were not on glucose-lowering medication. We carried forward CpGs identified in EWAS for T2D and meta-analysed EWAS for HbA1c and homeostatic model assessment estimates of insulin sensitivity (HOMA-S) as exposures to two-sample MR analysis. Independent cis methylation QTLs (meQTLs) were calculated using methylation data from blood and primary hepatocytes and subsequently used as instrumental variables (SNP-exposure associations). Genome-wide association analyses for T2D on 4,120 participants from the AADM study were used to derive the SNP-outcome associations. Longitudinal trait data (n=138) and RNAseq data (n=77 blood, 49 adipose, 55 muscle) available for a subset of Nigerians were used for follow-up analyses.

Results: We identified 28 CpGs associated with T2D, 26 with HbA1c, and 3 with HOMA-S (total 57), of which 49 had meQTLs in blood (AADM study data) and 4 had meQTLs in primary hepatocytes from African Americans. MR analysis provided evidence for causality for cg00036588 and cg16759041 in T2D using blood and hepatocyte meQTLs, respectively. Longitudinal analyses showed an association between baseline methylation of these CpGs with HbA1c at follow-up. RNAseq data revealed a cis correlation between cg00036588 with FAM83C (FDR = 3.3×10-4) and EIF6 (FDR = 0.13) in skeletal muscle.

Conclusions: Our study identified two epigenetic markers as likely causal for T2D in West African populations. In addition to enhancing our understanding of disease mechanisms, these CpGs with evidence of causal associations could be prioritized as potential biomarkers for early detection of disease or as drug development targets.




Come visit us at [CRGGH](https://www.genome.gov/about-nhgri/Center-for-Research-on-Genomics-and-Global-Health)!


## Structure <br>
The analysis is organized into three main phases:

---

### Phase 1. Epigenome-Wide Association Studies (EWAS)

Scripts in `1_EWAS/` include:

- `1a_cell_proportions.R`: Calculates estimated cell proportions using method by Houseman et al.  
- `1b_normalisation.R`: Filters probes, normalizes, and annotates to generate cleaned beta and M-value files 
- `1c_run_ewas.R`: Run EWAS on normalised M-values  
- `1d_meta-analysis.run`: Meta-analyse EWAS summary statistics from multiple cohorts using METAL 

---

### Phase 2. Mendelian Randomization (MR)

Scripts in `2_MR/` include:

- `2a_prepare_input.R`: Prepares Mvalue file for QTLtools
- `2b_calculate_meqtls.sh`: Calculates genome-wide meQTLs using QTLtools (SNP-exposure associations)
- `2c_prepare_instruments.sh`: Selects signficant meQTLs and performs LD clumping
- `2d_outcome_GWAS.R`: Runs GWAS on non-overlapping set of samples (SNP-outcome associations)
- `2e_run_mr_analysis.R`: Runs MR models including checks for MR assumptions

---

### Phase 3. Posthoc Analyses

Scripts in `3_Posthoc/` include:

- `3a_stability_cpgs.R`: Compares mean methylation betas between two timepoints  
- `3b_temporality_cpgs.R`: Plots baseline methylation by HbA1c at follow-up
- `3c_calculate_eqtms.R`: Calculates eQTMs using MatrixeQTL  

---

### Notes

- All file paths have been removed or replaced with relative paths.
- No filenames contain identifying information.

---

### Contact

For any questions, please contact:

**Karlijn A.C. Meeks** <br>
Assistant Professor
Division of Endocrinology, Diabetes and Nutrition
University of Maryland School of Medicine
670 W. Baltimore Street, HSF III 4061
Baltimore, MD 21201
E-mail: karlijn.meeks@som.umaryland.edu
Tel: (+1) 410 706 0132

Special Volunteer
Center for Research on Genomics and Global Health
National Human Genome Research Institute
National Institute of Health
E-mail: karlijn.meeks@nih.gov
