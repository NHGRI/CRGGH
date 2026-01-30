############################################################
# Title: 03b_METAL_AllWestAfricans.txt
# Project: Adiponectin EWAS – RODAM & AADM
# Author: Muhulo Muhau Mungamba
# Date: 2024-10-30
#
# Description:
# Fixed-effects meta-analysis of adiponectin EWAS results
# across RODAM and AADM, using stratified EWAS summary
# statistics (T2D cases and controls). This meta-analysis was performed using METAL and run
# in a Linux command-line environment (not in R)
#
# Input files were generated using the R script:
# 03a_prepare_input_files.R
############################################################

REMOVEFILTERS
SCHEME STDERR
USESTRAND OFF
GENOMICCONTROL OFF
AVERAGEFREQ OFF
MINMAXFREQ OFF
COLUMNCOUNTING LENIENT
SEPARATOR SPACE

# Column mappings
MARKER CPG
EFFECT BETA
PVALUE PVALUE
STDERR SE
WEIGHT N

# Output file
OUTFILE Adiponectin_AllWestAfricans_META

# Input files (prepared in 03a_prepare_input_files.R)
PROCESS RODAM_cases_controls_METAL.txt
PROCESS AADM_cases_controls_METAL.txt

# Run meta-analysis with heterogeneity statistics
ANALYZE HETEROGENEITY

######################################################
#END od analysis
#######################################################