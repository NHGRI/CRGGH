############################################################
# Title: 03f_METAL_AADM_all.R
# Project: Adiponectin EWAS – AADM
# Author: Muhulo Muhau Mungamba
# Date: 2024-10-30
#
# Description:
# Fixed-effects meta-analysis of adiponectin EWAS results
# within the AADM cohort (Nigerians), combining stratified
# EWAS summary statistics from T2D cases and T2D controls.
#
# This meta-analysis was performed using METAL and run
# in a Linux command-line environment (not in R).
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

OUTFILE Adiponectin_AADM_All_META

# -------------------------------
# Input files (AADM only:
# T2D cases + T2D controls)
# -------------------------------

PROCESS AADM_T2Dcases_METAL.txt
PROCESS AADM_T2Dcontrols_METAL.txt

# -------------------------------
# Run meta-analysis with
# heterogeneity statistics
# -------------------------------

ANALYZE HETEROGENEITY

############################################################
# END OF SCRIPT
############################################################
