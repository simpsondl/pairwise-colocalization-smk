library(dplyr)

# Logging helper
log_msg <- function(...) {
  msg <- paste(..., collapse = "")
  timestamped_msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  message(timestamped_msg)
  if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
    cat(timestamped_msg, file = snakemake@log[[1]], append = TRUE, sep = "\n")
  }
}

# Read pairs file
log_msg("Reading pairs file: ", snakemake@input[["pairs"]])
pairs <- read.csv(snakemake@input[["pairs"]])
log_msg("Loaded ", nrow(pairs), " pairs from input file")

# Get info for this pair
log_msg("Looking for pair_id: ", snakemake@params[["pair_id"]])
pair <- pairs[pairs$pair_index == snakemake@params[["pair_id"]], ]

if(nrow(pair) == 0) {
  log_msg("ERROR: Pair ID not found in pairs file: ", snakemake@params[["pair_id"]])
  stop("Pair ID not found in pairs file: ", snakemake@params[["pair_id"]])
}

log_msg("Found pair: ", pair$gwas1, " vs ", pair$gwas2, " on chr", pair$chr, " (", pair$pos.start, "-", pair$pos.end, ")")
log_msg("Starting extract_regions.R")

# Function to check if file has a header
check_for_header <- function(file_path, expected_cols = 7) {
  # Read first few lines to inspect
  first_lines <- readLines(file_path, n = 5)
  
  if(length(first_lines) == 0) {
    log_msg("WARNING: File is empty: ", file_path)
    return(list(has_header = FALSE, first_row = character(0)))
  }
  
  # Split first line by tabs/spaces
  first_row <- strsplit(first_lines[1], "\t| +")[[1]]
  
  # Check if first row has expected number of columns
  if(length(first_row) != expected_cols) {
    log_msg("WARNING: First row has ", length(first_row), " columns, expected ", expected_cols)
  }
  
  # Check if first row looks like header (contains common GWAS column names)
  header_keywords <- c("chr", "chrom", "pos", "position", "snp", "rs", "allele", "beta", "se", "pval", "p.value", "p_val")
  has_header_keywords <- any(sapply(header_keywords, function(keyword) {
    any(grepl(keyword, first_row, ignore.case = TRUE))
  }))
  
  # Check if first row looks numeric (suggests no header)
  is_numeric_row <- all(sapply(first_row, function(x) {
    !is.na(suppressWarnings(as.numeric(x)))
  }))
  
  has_header <- has_header_keywords && !is_numeric_row
  
  log_msg("Header analysis for ", basename(file_path), ": has_header=", has_header, 
          ", first_row_looks_numeric=", is_numeric_row, ", contains_header_keywords=", has_header_keywords)
  
  return(list(has_header = has_header, first_row = first_row))
}

# Function to extract region from a GWAS file
extract_region <- function(gwas_file, chr, start, end) {
  log_msg("Extracting region chr", chr, ":", start, "-", end, " from ", basename(gwas_file))
  
  if(!file.exists(gwas_file)){
    log_msg("WARNING: GWAS file not found: ", gwas_file)
    return(data.frame())
  }
  
  # Check for header
  header_check <- check_for_header(gwas_file)
  use_header <- header_check$has_header
  
  log_msg("Reading GWAS file with header=", use_header)
  data <- read.table(gwas_file, header = use_header)
  
  if(use_header) {
    log_msg("File has header, columns: ", paste(colnames(data), collapse = ", "))
  } else {
    log_msg("File has no header, assigning standard column names")
    colnames(data) <- c("chromosome",	"position",	"allele1",	"allele2",	"beta",	"se",	"pval")
  }
  
  log_msg("Read ", nrow(data), " rows from GWAS file")
  
  # Extract region
  log_msg("Filtering for chromosome ", chr, " and position range ", start, "-", end)
  region_data <- data[data$chromosome == chr & 
                     data$position >= start & 
                     data$position <= end, ]
  
  log_msg("Extracted ", nrow(region_data), " SNPs in the target region")
  return(region_data)
}

# Extract regions for both GWAS
log_msg("Constructing file paths for GWAS summary statistics")
gwas1_file <- file.path(snakemake@input[["sumstats_dir"]], paste0(pair$gwas1, ".txt"))
gwas2_file <- file.path(snakemake@input[["sumstats_dir"]], paste0(pair$gwas2, ".txt"))

log_msg("Extracting region for GWAS1: ", pair$gwas1)
gwas1 <- extract_region(gwas1_file, pair$chr, pair$pos.start, pair$pos.end)
log_msg("Extracting region for GWAS2: ", pair$gwas2)
gwas2 <- extract_region(gwas2_file, pair$chr, pair$pos.start, pair$pos.end)

if(nrow(gwas1) == 0 || nrow(gwas2) == 0){
  log_msg("WARNING: One or both GWAS region files are empty")
  log_msg("GWAS1 rows: ", nrow(gwas1), ", GWAS2 rows: ", nrow(gwas2))
  log_msg("One or both GWAS region files are empty; writing empty output: ", snakemake@output[["output_joined_sumstats"]])
  write.csv(data.frame(), snakemake@output[["output_joined_sumstats"]], row.names=FALSE, quote=FALSE)
} else {
  log_msg("Performing inner join on chromosome, position, and alleles")
  both_gwas <- inner_join(gwas1[,1:6], gwas2[,1:6], 
                          by=c("chromosome", "position", "allele1", "allele2"), 
                          suffix=c(".gwas1", ".gwas2"))
  log_msg("Joined dataset contains ", nrow(both_gwas), " rows")
  log_msg("Writing joined sumstats to: ", snakemake@output[["output_joined_sumstats"]])
  write.csv(both_gwas, snakemake@output[["output_joined_sumstats"]], row.names=FALSE, quote=FALSE)
}

log_msg("extract_regions.R finished for pair_id=", pair$pair_index)