############################################################
# Title: 03e_METAL_RODAM_all.R
# Project: Adiponectin EWAS – RODAM
# Author: Muhulo Muhau Mungamba
# Date: 2024-10-30
#
# Description:
# Fixed-effects meta-analysis of adiponectin EWAS results
# within the RODAM cohort (Ghanaians), combining stratified
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

OUTFILE Adiponectin_RODAM_All_META

# -------------------------------
# Input files (RODAM only:
# T2D cases + T2D controls)
# -------------------------------

PROCESS RODAM_T2Dcases_METAL.txt
PROCESS RODAM_T2Dcontrols_METAL.txt

# -------------------------------
# Run meta-analysis with
# heterogeneity statistics
# -------------------------------

ANALYZE HETEROGENEITY

############################################################
# END OF SCRIPT
############################################################
