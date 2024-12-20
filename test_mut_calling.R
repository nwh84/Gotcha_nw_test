
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
source("/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_nw_test/R/FastqFilteringNew.R")

print("start time:")
print(Sys.time())

path <- "/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_test_pipeline/GoTChA_test_data2/aa"

out <- MutationCalling(out = path,
                    which.read = "R1", # Look into R1 now
                    mutation.start = 1, # take the whole read
                    mutation.end = 50, # take the whole read
                    barcodes.file.path = "/sc/arion/projects/MDS/noelle/gotcha_pipeline_nwheeler/Gotcha_test_pipeline/GoTChA_test_data2/singlecell.csv", # barcodes found in the experiment
                    primer.sequence = "A", # this checks the match for the intial nucleotides
                    primed.max.mismatch = 1, # putting 1 here allows one nucleotide mismatch a thus primer matching should not eliminate any reads as we only put 1 nucleotide anyway
                    wt.sequence = "TACAGTGCAGGGGAAAGAATAGTAGACATAATAGCAACAGACATACAAAC", # Pol
                    mut.sequence= "CATCTGGCCTGGTGCAATAGGCCCTGCATGCACTGGATGCAATCTATCCC", # Pol as well, second most abundant
                    wt.max.mismatch = 0, # set up on perfect matching for increased accuracy
                    mut.max.mismatch = 0, # set up on perfect matching for increased accuracy
                    ncores = 1,
                    testing = FALSE)




print("finished time:")
print(Sys.time())

saveRDS(out, file = "output_test2")


