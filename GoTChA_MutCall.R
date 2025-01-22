
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
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/MutationCallingNew.R")
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/FastqFiltering.R")

start.time <- Sys.time()

#path <- "/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_test_pipeline/GoTChA_test_data/"

path <- "/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_test_pipeline/GoTChA_test_data/gotcha_test_data_subset2/"

BatchMutationCalling(out  = path,
                     which.read = "R1",
                     mutation.start = 1, # take the whole read
                     mutation.end = 50,
                     barcodes.file.path = "/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_test_pipeline/GoTChA_test_data/gotcha_test_data_subset/CBNN_0_singlecell.csv",
                     wt.sequence = "GCATGTATGCAATGCCTTGGTAGGAATGGGACAGGTGTAGGATGGAAAAT",
                     mut.sequence= "AAGGCGTTTCTTCTCTGACCGCACAACTGGGGCCTGGGGGGCTCCAAAGC",
                     wt.max.mismatch = 0, # set up on perfect matching for increased accuracy
                     mut.max.mismatch = 0,
                     ncores = 44,
                     rewrite = FALSE,
                     testing = TRUE,
                     max.distance = 3
)


batch.mutation.time <- Sys.time()

MergeMutationCalling(out = path)

end.time <- Sys.time()

print("start time:")
print(start.time)
print("finished batch mutation:")
print(batch.mutation.time)
print("finished merging files:")
print(end.time)


