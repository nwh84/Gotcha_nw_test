#' Define read genotype and read counts per genotype for each cell barcode
#' @param out Path to the fastq or filtered fastq files
#' @param wt.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the wild-type sequence
#' @param mut.max.mismatch Integer indicating the number of accepted missmatches when performing pattern matching for the mutant sequence
#' @param keep.raw.reads Logical. Whether to return the raw reads in the output file. Defaults to false
#' @param ncores Integer indicating the number of cores to use for parallel processing
#' @param reverse.complement Whether to take the reverse complement of the cell barcodes
#' @param testing Logical indicating whether to sample the first 1,000 reads for testing the function
#' @param which.read Which read to select to look for the mutation site
#' @param primer.sequence Character vector of length one indicating the primer sequence
#' @param primed.max.mismatch  Integer indicating the maximum number of mismatches accepted when searching for the primer sequence
#' @param barcodes.file.path Path to the file containing the cell barcodes detected in the experiment
#' @param wt.sequence Character vector of length one specifying the expected wild-type sequence
#' @param mut.sequence Character vector of length one specifying the expected mutant sequence
#' @param mutation.start Position in which the expected wild-type or mutant sequence starts in the read
#' @param mutation.end Position in which the expected wild-type or mutant sequence ends in the read
#' @param max.distance Maximum number of mismatches allowed between barcodes and whitelist
#' @return output Datatable with barcode and genotype calls


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

convert_to_numeric_matrix_raw <- function(input) {
  # Convert each character in the barcode to its ASCII value
  as.data.table(do.call(rbind, lapply(input, function(bc) {
    as.integer(charToRaw(bc))
  })))
}

convert_to_numeric_matrix <- function(input, num) {
  char_to_num <- c(A = 1, C = 2, G = 3, T = 4)
  numeric_matrix <- t(sapply(input, function(x) {
         char_to_num[unlist(strsplit(x, ""))]
     }))
  dimnames(numeric_matrix) <- NULL
  numeric_matrix <- as.data.table(numeric_matrix)
  return(numeric_matrix)
}

convert_to_numeric_matrix_reorder <- function(input, num) {
  # Convert each character in the barcode to its ASCII value
  #as.data.table(do.call(rbind, lapply(input, function(bc) {
  #  as.integer(charToRaw(bc))
  #})))
  if (num == 1){
    char_to_num <- c(A = 1, T = 2, C = 3, G = 4)
  } else if(num ==2){
    char_to_num <- c(T = 1, A = 2, G = 3, C = 4)
  } else if(num ==3){
    char_to_num <- c(C = 1, G = 2, A = 3, T = 4)
  } else {
    char_to_num <- c(T = 1, C = 2, A = 3, G = 4)
  }
  
  numeric_matrix <- t(sapply(input, function(x) {
    char_to_num[unlist(strsplit(x, ""))]
  }))
  dimnames(numeric_matrix) <- NULL
  numeric_matrix <- as.data.table(numeric_matrix)
  return(numeric_matrix)
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

# input is list of barcodes, whitelist and number of nearest neighbor desired
# output a matrix where each row 
rann_matching <- function(barcodes, whitelist, order, radius, nearest_neighbors) {
  # Convert barcodes and whitelist to numeric matrices
  if (order ==1){
    barcode_matrix <- convert_to_numeric_matrix(barcodes)
    whitelist_matrix <- convert_to_numeric_matrix(whitelist)
  } else{
    barcode_matrix <- convert_to_numeric_matrix_reorder(barcodes)
    whitelist_matrix <- convert_to_numeric_matrix_reorder(whitelist)
  }
  # Use RANN for nearest neighbor search
  nn_results <- nn2(
    data = whitelist_matrix,    # Whitelist data
    query = barcode_matrix,     # Query barcodes
    k = nearest_neighbors,      # Find num of nearest neighbor
    searchtype = "radius",      # Use radius-based search
    radius = radius       # Maximum allowed distance, will only allow for 2 mismatches
  )
  whitelist_ind <- nn_results$nn.idx
  # convert 0 to NA
  whitelist_ind <- replace(whitelist_ind, whitelist_ind==0, NA)
  # return the nearest whitelist barcodes
  mapped_matrix <- matrix(whitelist[whitelist_ind], nrow = nrow(nn_results$nn.idx), ncol = ncol(nn_results$nn.idx))
  return(mapped_matrix)
}

# for each barcode will output closest match to whitelist, "no match", or NA if there are two matches
# input is a barcodes, a list a whitelist barcodes, and max.distance allowed
hamming_matching <- function(barcodes, whitelist, max.distance) {
  # get number of mismatches between strings
  dist_mat <- stringdistmatrix(barcodes, whitelist, method = "hamming")
  num_mismatch <- apply(dist_mat, 1, min) # will find smallest distance, NAs are already removed
  # check if all whitelist are NA
  if(is.na(num_mismatch)){
    return("No Match")
    } else if (num_mismatch > max.distance) {
      return("No Match") # too many mismatches 
      } else if(length(which(dist_mat<=num_mismatch)) > 1) {
        return("No Match") # ambiguous mapping, should never happend
        } else {
          # take item in whitelist with lowest number of mismatch
          return(whitelist[which(dist_mat<=num_mismatch)])
  }
  return("No Match")
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
  if (testing) fastq_data <- subset_for_testing(fastq_data, max_reads = 10000, ncores = ncores)
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
  
  # print percent barcodes left, and run time for each 
  # run RANN to get 2 nearest neighbors
  non_match_barcodes <- as.character(barcodes$x)[!matched_barcodes_ind]
  matched_barcodes_rann <- rann_matching(barcodes=non_match_barcodes, whitelist = whitelist, order = 1, radius = 2.3, nearest_neighbors = 2)
  # run hamming apply on each row of matched barcodes
  result <- unlist(lapply(seq_along(non_match_barcodes), function(i) {
    hamming_matching(non_match_barcodes[i], matched_barcodes_rann[i], max.distance)
  }))
  # add the imperfect match barcodes
  matched_barcodes[!matched_barcodes_ind] <- result
  
  # run rann with reorder = 2
  ####
  final_non_match_barcodes <- as.character(barcodes$x)[which(matched_barcodes == "No Match")]
  final_matched_barcodes_rann <- rann_matching(barcodes=final_non_match_barcodes, whitelist = whitelist, order = 2, radius = 2.3, nearest_neighbors = 2)
  result <- unlist(lapply(seq_along(final_non_match_barcodes), function(i) {
    hamming_matching(final_non_match_barcodes[i], final_matched_barcodes_rann[i], max.distance)
  }))
  matched_barcodes[which(matched_barcodes == "No Match")] <- result

  # run rann with reorder = 3
  ####
  final_non_match_barcodes <- as.character(barcodes$x)[which(matched_barcodes == "No Match")]
  final_matched_barcodes_rann <- rann_matching(barcodes=final_non_match_barcodes, whitelist = whitelist, order = 3, radius = 2.3, nearest_neighbors = 2)
  result <- unlist(lapply(seq_along(final_non_match_barcodes), function(i) {
    hamming_matching(final_non_match_barcodes[i], final_matched_barcodes_rann[i], max.distance)
  }))
  matched_barcodes[which(matched_barcodes == "No Match")] <- result

# run rann with reorder = 4
  ####
  final_non_match_barcodes <- as.character(barcodes$x)[which(matched_barcodes == "No Match")]
  final_matched_barcodes_rann <- rann_matching(barcodes=final_non_match_barcodes, whitelist = whitelist, order = 4, radius = 2.3, nearest_neighbors = 2)
  result <- unlist(lapply(seq_along(final_non_match_barcodes), function(i) {
    hamming_matching(final_non_match_barcodes[i], final_matched_barcodes_rann[i], max.distance)
  }))
  matched_barcodes[which(matched_barcodes == "No Match")] <- result
  ####
  
  # run hamming on the remaining no match barcodes
  #final_non_match_barcodes <- as.character(barcodes$x)[which(matched_barcodes == "No Match")]
  #result_last_barcodes <- unlist(lapply(seq_along(final_non_match_barcodes), function(i) {
  #  hamming_matching(final_non_match_barcodes[i], whitelist, max.distance)
  #}))
  # add the final hamming barcodes for those that didn't have match
  #matched_barcodes[which(matched_barcodes == "No Match")] <- result_last_barcodes
  
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
