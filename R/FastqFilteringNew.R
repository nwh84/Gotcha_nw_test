FastqFiltering <- function(
    path, 
    out = "/path_to_fastqs/", 
    reads = 2000000, 
    min.quality = 15, 
    min.bases = 1, 
    which.read = "R1", 
    read.region = NULL, 
    ncores = 1
) {
  library(ShortRead)
  library(parallel)
  
  message("------- STARTING FASTQSPLIT -------")
  
  # Create output directories if they don't exist
  split_dir <- file.path(out, "Split")
  filtered_dir <- file.path(out, "Filtered")
  dir.create(split_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(filtered_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Find input FASTQ files
  file.index <- list.files(path, pattern = "\\.fastq\\.gz$", full.names = TRUE)
  if (length(file.index) == 0) stop("No fastq files detected in path folder")
  
  message("------- SPLITTING FASTQS -------")
  mclapply(file.index, function(x) {
    cmd <- sprintf(
      "gunzip -c %s | split - -l %d --filter='gzip -f > $FILE.gz' %s/Split_$(basename %s .fastq.gz).fastq",
      x, reads * 4, split_dir, x
    )
    system(cmd)
  }, mc.cores = ncores)
  
  message("------- FASTQS SPLIT COMPLETED -------")
  
  # Get split FASTQ files
  fastq.files <- list.files(split_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE)
  if (length(fastq.files) == 0) stop("No split FASTQ files detected")
  
  # Filter files based on read type
  selected.files <- grep(which.read, fastq.files, value = TRUE)
  if (length(selected.files) == 0) stop(sprintf("No FASTQ files matched '%s'", which.read))
  
  message("------- STARTING FILTERING -------")
  
  # Define filtering function
  filter_fastq <- function(file) {
    fq <- readFastq(file)
    qual_matrix <- as(quality(fq), "matrix")
    
    # Apply region filter if specified
    if (!is.null(read.region)) {
      qual_matrix <- qual_matrix[, read.region, drop = FALSE]
    }
    
    # Apply quality and length filtering
    valid_reads <- rowSums(qual_matrix <= min.quality) < min.bases
    fq[valid_reads]
  }
  
  # Filter each FASTQ file
  mclapply(selected.files, function(file) {
    filtered <- filter_fastq(file)
    output_file <- file.path(filtered_dir, paste0("filtered_", basename(file)))
    writeFastq(filtered, output_file)
  }, mc.cores = ncores)
  
  message("------- FILTERING COMPLETED -------")
  
  # Remove temporary split files
  message("------- CLEANING UP TEMPORARY FILES -------")
  unlink(split_dir, recursive = TRUE)
  message("Temporary split files removed.")
  
  # Metrics
  filtered_files <- list.files(filtered_dir, pattern = "\\.fastq$", full.names = TRUE)
  total_reads <- unlist(mclapply(filtered_files, function(file) {
    length(readFastq(file))
  }, mc.cores = ncores))
  
  message(sprintf("Mean number of filtered reads: %.2f", mean(total_reads)))
  message(sprintf("Percentage of remaining reads: %.2f%%", 
                  mean(total_reads) / (reads * length(file.index)) * 100))
  
  message("------- FASTQ FILTERING COMPLETE -------")
}