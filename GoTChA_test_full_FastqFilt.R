
.libPaths(c("/hpc/packages/minerva-rocky9/rpackages/4.2.0/site-library", "/hpc/packages/minerva-rocky9/rpackages/bioconductor/3.15", "/hpc/packages/mixnerva-rocky9/R/4.2.0/1ib64/R/library"))

library(dplyr)
library(parallel)
library(readr)
library(reticulate)
library(lme4)
library(methods)
library(stats)
library(utils)
library(stringdist)
library(ShortRead)
library(Biostrings)
library(graphics)
library(RANN)
library(data.table)

source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/AddGenotypingArchr.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/BatchMutationCalling.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/DiffLMM.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/DiffPeaks.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/FilterGenotyping.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/GotchaLabeling.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/MergeFastqs.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/MergeMutationCalling.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/MutationCalling.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/FastqFiltering.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/FastqSplit.R")


path <- "/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_test_pipeline/GoTChA_test_data/"

FastqFiltering(out = path, min.quality = 15, which.read = "R1", ncores = 1, read.region = 15, min.bases= 5)


