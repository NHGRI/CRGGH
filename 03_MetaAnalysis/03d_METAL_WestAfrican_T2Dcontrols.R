############################################################
# Title: 03d_METAL_WestAfrican_T2Dcontrols.txt
# Project: Adiponectin EWAS – RODAM & AADM
# Author: Muhulo Muhau Mungamba
# Date: 2024-10-30
#
# Description:
# Fixed-effects meta-analysis of adiponectin EWAS results
# among participants without type 2 diabetes (T2D controls),
# combining stratified EWAS summary statistics from the
# RODAM (Ghanaians) and AADM (Nigerians) cohorts.
#
# This meta-analysis was performed using METAL and run
# in a Linux command-line environment (not in R)
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

OUTFILE Adiponectin_WestAfricans_T2Dcontrols_META

# -------------------------------
# Input files (T2D controls only)
# -------------------------------

PROCESS RODAM_T2Dcontrols_METAL.txt
PROCESS AADM_T2Dcontrols_METAL.txt

# -------------------------------
# Run meta-analysis with
# heterogeneity statistics
# -------------------------------

ANALYZE HETEROGENEITY

############################################################
# END OF SCRIPT
############################################################
