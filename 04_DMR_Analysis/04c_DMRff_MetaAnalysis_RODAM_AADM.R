############################################
# Title: Meta-analysis of Differentially Methylated Regions (DMRs)
# Method: DMRff
# Cohorts: RODAM (Ghanaians) and AADM (Nigerians)
# Outcome: Adiponectin EWAS
#
# Author: Muhulo Muhau Mungamba
# Date: 2025-10-17
#
# Description:
# This script performs region-level meta-analysis of DNA methylation
# using EWAS summary statistics from RODAM and AADM.
#
# DMRs are identified using DMRff, which leverages spatial correlation
# between adjacent CpGs and combines regression coefficients and
# standard errors across cohorts.
#
# This analysis corresponds to the DMRff meta-analysis described
# in the Methods section of the manuscript.
############################################

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)

# Load required libraries
library(dmrff)
library(dplyr)
library(GenomicRanges)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

### ------------------------------------------------------------------
### Load EWAS summary statistics
### ------------------------------------------------------------------

# RODAM summary statistics (combined T2D cases + controls)
rodam <- read.delim(
  "/project/MUHAU_RODAM/EWAS_Adiponectin/Results/TopTables/2025.10.14_TopTable_RODAM_Adiponectin_combined_T2Dcasescontrols_n315_age_sex_cells_pos_batcharry_BMI_site_DM.txt",
  sep = "\t",
  header = TRUE
)

# AADM summary statistics (combined T2D cases + controls)
aadm <- read.delim(
  "/project/MUHAU_RODAM/EWAS_Adiponectin/Results/TopTables/2024.10.22_TopTable_AADM_Adiponectin_n593_age_sex_cells_bmi_t2d.txt",
  sep = "\t",
  header = TRUE
)

### ------------------------------------------------------------------
### Harmonize columns for DMRff
### ------------------------------------------------------------------

rodam_ready <- rodam %>%
  transmute(
    CpG  = CpG,
    chr  = chr,
    pos  = pos,
    beta = logFCBeta,
    se   = SE_mval,
    pval = P.Value.Mval
  ) %>%
  filter(!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se))

aadm_ready <- aadm %>%
  transmute(
    CpG  = CpG,
    chr  = chr,
    pos  = pos,
    beta = logFCBeta,
    se   = SE_beta,
    pval = P.Value
  ) %>%
  filter(!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se))

### ------------------------------------------------------------------
### Restrict to CpGs shared across cohorts
### ------------------------------------------------------------------

shared_cpgs <- intersect(rodam_ready$CpG, aadm_ready$CpG)

rodam_ready <- rodam_ready %>% filter(CpG %in% shared_cpgs)
aadm_ready  <- aadm_ready  %>% filter(CpG %in% shared_cpgs)

### ------------------------------------------------------------------
### Prepare DMRff inputs
### ------------------------------------------------------------------

pre_rodam <- dmrff.pre(
  estimate = rodam_ready$beta,
  se       = rodam_ready$se,
  chr      = rodam_ready$chr,
  pos      = rodam_ready$pos
)

pre_aadm <- dmrff.pre(
  estimate = aadm_ready$beta,
  se       = aadm_ready$se,
  chr      = aadm_ready$chr,
  pos      = aadm_ready$pos
)

meta_input <- list(pre_rodam, pre_aadm)

### ------------------------------------------------------------------
### Run DMRff meta-analysis
### ------------------------------------------------------------------

dmr_meta <- dmrff.meta(meta_input)

dmrs <- dmr_meta$dmrs

# Keep statistically significant DMRs
sig_dmrs <- dmrs %>%
  filter(p.adjust < 0.05 & n >= 2) %>%
  arrange(p.adjust)

### ------------------------------------------------------------------
### Annotate DMRs with genes
### ------------------------------------------------------------------

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- as.data.frame(anno)[, c("Name", "chr", "pos", "UCSC_RefGene_Name")]

dmr_ranges <- GRanges(
  seqnames = sig_dmrs$chr,
  ranges   = IRanges(start = sig_dmrs$start, end = sig_dmrs$end)
)

anno_ranges <- GRanges(
  seqnames = anno_df$chr,
  ranges   = IRanges(start = anno_df$pos, end = anno_df$pos),
  gene     = anno_df$UCSC_RefGene_Name
)

hits <- findOverlaps(dmr_ranges, anno_ranges)

annotated_dmrs <- sig_dmrs[queryHits(hits), ]
annotated_dmrs$Gene <- anno_df$UCSC_RefGene_Name[subjectHits(hits)]

### ------------------------------------------------------------------
### Save outputs
### ------------------------------------------------------------------

write.csv(
  annotated_dmrs,
  "2025.03.25_DMRff_MetaAnalysis_RODAM_AADM_significant_DMRs.csv",
  row.names = FALSE
)

############################################
# END OF SCRIPT
############################################
