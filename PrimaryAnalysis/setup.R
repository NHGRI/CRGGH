# Setup notes -------------------------------------------------------------

# setup.R referenced by other scripts to avoid having to copy/paste/maintain this section for every analysis

# Global parameters -------------------------------------------------------

rnaProject <- "skmus"
regression.vars <- c("SampleSex", "SampleAge", "sequencerID", "nCount_RNA")
cum.var.thresh <- 90
resolution <- 0.5
do.sctransform <- "each" # one of FALSE, each, pooled

## infrequently modified
do.doubletFinder <- TRUE
min.cells <- 3
min.features <- 150
doublet.var.thresh <- 90
predicted.doubletRate <- 0.05


# Directories -------------------------------------------------------------

rna.dir <- "."
path_to_data <- "./preprocessing"
sourceable.functions <- list.files(path = "./RFunctions", pattern = "*.R", full.names = TRUE)
metadata <- "./MetaData.txt"

# Load libraries ----------------------------------------------------------

library(dplyr)
library(Seurat)
library(patchwork)
library(cowplot)
library(foreach)
library(doParallel)

# load local functions
sourceable.functions <- sourceable.functions[grep(c("ColorPalette.R|PercentVariance.R|runDoubletFinder.R|AssignMetadata.R"), sourceable.functions)]
invisible(sapply(sourceable.functions, source))
color.palette <- c(color.palette, color.palette)

# create output files
try(setwd(rna.dir), silent = TRUE)

# define parallel
if(!grepl("Darwin", comp.type, ignore.case = TRUE)){
  options(future.globals.maxSize = 200000 * 1024^2)
  cl <- makeCluster(future::availableCores(), outfile = "")
  registerDoParallel(cl)
}

##session info
session.file <- paste0(rnaProject, "-analysis_sessionInfo.txt")
write(Sys.time(),file=session.file,append = TRUE)
write(capture.output(sessionInfo()), session.file, append = TRUE)
