library(data.table)
library(parallel)
library(ShortRead)
library(RANN)

# Helper functions
read_and_process_fastq <- function(path, pattern, ncores) {
  fastq_files <- dir(path, pattern = pattern, full.names = TRUE)
  if (length(fastq_files) == 0) stop("No fastq files detected.")
  res = mclapply(fastq_files, function(file) {
    temp <- readFastq(file)
    as.data.table(sread(temp))
  }, mc.cores = ncores)
  names(res) = fastq_files
  return(res)
}

subset_for_testing <- function(fastq_data, max_reads = 1000, ncores) {
  mclapply(fastq_data, function(dt) {
    if (nrow(dt) > max_reads) dt[1:max_reads, ] else dt
  }, mc.cores = ncores)
}

convert_to_numeric_matrix <- function(input) {
  # Convert each character in the barcode to its ASCII value
  as.data.table(do.call(rbind, lapply(input, function(bc) {
    as.integer(charToRaw(bc))
  })))
}

load_whitelist <- function(file_path) {
  if (!file.exists(file_path)) stop("Whitelist file not found: ", file_path)
  whitelist_data <- fread(file_path)
  if ("barcode" %in% colnames(whitelist_data)) {
    whitelist <- whitelist_data$barcode
  } else {
    # Assume a single-column CSV without a header
    whitelist <- whitelist_data[[1]]
  }
  # Trim to consistent length if needed (e.g., 16 characters for 10x barcodes)
  whitelist <- substr(whitelist, 1, 16)
  return(unique(whitelist))
}

rann_matching <- function(barcodes, whitelist, max.distance) {
  # Convert barcodes and whitelist to numeric matrices
  barcode_matrix <- convert_to_numeric_matrix(barcodes)
  whitelist_matrix <- convert_to_numeric_matrix(whitelist)
  # Use RANN for nearest neighbor search
  nn_results <- nn2(
    data = whitelist_matrix,    # Whitelist data
    query = barcode_matrix,     # Query barcodes
    k = 1,                      # Find 1 nearest neighbor
    searchtype = "radius",      # Use radius-based search
    radius = max.distance       # Maximum allowed distance
  )
  # Extract results
  matched_indices <- nn_results$nn.idx
  distances <- nn_results$nn.dists
  # Map matches back to whitelist or provide no-match labels
  matches <- sapply(seq_along(barcodes), function(i) {
    if (matched_indices[i] == 0 || distances[i] > max.distance) {
      return("No Match")
    } else {
      return(whitelist[matched_indices[i]])
    }
  })
  return(matches)
}

genotype_reads <- function(reads, wt_seq, mut_seq, mutation_start, mutation_end, wt_max_mismatch, mut_max_mismatch, ncores) {
  mclapply(reads, function(read) {
    wt_count <- vcountPattern(wt_seq, substr(read, mutation_start, mutation_end), max.mismatch = wt_max_mismatch)
    mut_count <- vcountPattern(mut_seq, substr(read, mutation_start, mutation_end), max.mismatch = mut_max_mismatch)
    genotype <- ifelse(wt_count == 1 & mut_count == 1, "Ambiguous",
                       ifelse(wt_count == 0 & mut_count == 0, "No information",
                              ifelse(wt_count == 1 & mut_count == 0, "WT", "MUTANT")))
    data.table(WT = wt_count, MUT = mut_count, Genotype = genotype)
  }, mc.cores = ncores)
}
# Main function
MutationCalling <- function(out, barcodes.file.path, wt.max.mismatch = 0, mut.max.mismatch = 0,
                            keep.raw.reads = FALSE, ncores = 1, reverse.complement = TRUE,
                            testing = FALSE, which.read = "R1", primer.sequence = "CCTCATCATCCTCCTTGTC",
                            primed.max.mismatch = 3, wt.sequence = "CGG", mut.sequence = "CAG",
                            mutation.start = 31, mutation.end = 34, max.distance = 2) {
  message("------- BEGIN MUTATION CALLING -------")
  # Load whitelist from the specified file path
  whitelist <- load_whitelist(barcodes.file.path)
  # remove "NO_BARCODE" from whitelist
  whitelist <- whitelist[whitelist != "NO_BARCODE"]
  message("Whitelist loaded with ", length(whitelist), " unique barcodes.")
  # Load FASTQ files
  fastq_data <- read_and_process_fastq(out, pattern = ".fastq.gz", ncores = ncores)
  if (testing) fastq_data <- subset_for_testing(fastq_data, max_reads = 1000, ncores = ncores)
  message("------- FASTQ FILES LOADED -------")
  # Process sequences
  barcodes <- fastq_data[[grep(names(fastq_data), pattern = "_R2_")]]
  reads <- fastq_data[[grep(names(fastq_data), pattern = paste0("_", which.read, "_"))]]
  # reverse complement and convert to dnastringset
  if (reverse.complement) {
    barcodes <- lapply(barcodes, function(x) reverseComplement(DNAStringSet(x)))
    message("Cell barcodes have been reverse complemented.")
  }

  message("------- STARTING BARCODE MATCHING -------")
  # check which barcodes are perfect match
  matched_barcodes_ind <- (as.character(barcodes$x) %in% whitelist)
  matched_barcodes <- rep(NA, length(matched_barcodes_ind))
  # add the perfect match barcodes
  matched_barcodes[matched_barcodes_ind] <- as.character(barcodes$x)[matched_barcodes_ind]

  # Barcode matching against the whitelist using RANN
  matched_barcodes_rann <- rann_matching(as.character(barcodes$x)[!matched_barcodes_ind], whitelist, max.distance)
  # add the imperfect match barcodes
  matched_barcodes[!matched_barcodes_ind] <- matched_barcodes_rann

  message("------- BARCODE MATCHING COMPLETED -------")
  message("------- STARTING PER READ GENOTYPING -------")
  # Genotype reads
  genotyped_reads <- genotype_reads(reads, wt.sequence, mut.sequence, mutation.start, mutation.end, wt.max.mismatch, mut.max.mismatch, ncores)
  message("------- GENOTYPING COMPLETED -------")
  message("------- SAVING OUTPUT... -------")
  # Output processing
  output <- list(matched_barcodes = matched_barcodes, genotyped_reads = genotyped_reads)
  if (keep.raw.reads) output$raw_reads <- reads
  message("DONE!")
  return(output)
}
