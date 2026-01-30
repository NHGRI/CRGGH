############################################################
# Title: 03c_METAL_WestAfrican_T2Dcases.txt
# Project: Adiponectin EWAS – RODAM & AADM
# Author: Muhulo Muhau Mungamba
# Date: 2024-10-30
#
# Description:
# Fixed-effects meta-analysis of adiponectin EWAS results
# among type 2 diabetes (T2D) cases only, combining
# stratified EWAS summary statistics from the RODAM
# (Ghanaians) and AADM (Nigerians) cohorts.
#
# This meta-analysis was performed using METAL and run
# in a Linux command-line environment (not in R), exactly
# as described in the eBioMedicine manuscript.
#
# Input files were generated using the R script:
# 03a_prepare_input_files.R
############################################################

# -------------------------------
# METAL configuration
# -------------------------------

REMOVEFILTERS
SCHEME STDERR
USESTRAND OFF
GENOMICCONTROL OFF
AVERAGEFREQ OFF
MINMAXFREQ OFF
COLUMNCOUNTING LENIENT
SEPARATOR SPACE

# -------------------------------
# Column mappings
# -------------------------------

MARKER CPG
EFFECT BETA
PVALUE PVALUE
STDERR SE
WEIGHT N

# -------------------------------
# Output file
# -------------------------------

OUTFILE Adiponectin_WestAfricans_T2Dcases_META

# -------------------------------
# Input files (T2D cases only)
# -------------------------------

PROCESS RODAM_T2Dcases_METAL.txt
PROCESS AADM_T2Dcases_METAL.txt

# -------------------------------
# Run meta-analysis with
# heterogeneity statistics
# -------------------------------

ANALYZE HETEROGENEITY

############################################################
# END OF SCRIPT
############################################################
