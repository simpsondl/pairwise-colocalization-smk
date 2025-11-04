# Read input signals
log_msg <- function(...) {
  msg <- paste(..., collapse = "")
  timestamped_msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  message(timestamped_msg)
  if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
    cat(timestamped_msg, file = snakemake@log[[1]], append = TRUE, sep = "\n")
  }
}

log_msg("Reading input signals file")
signals <- read.csv(snakemake@input[["signals"]])
log_msg(sprintf("Loaded %d signals from input file", nrow(signals)))

pairs <- data.frame()

# Find pairs within distance threshold
distance <- snakemake@params[["distance"]]
log_msg(sprintf("Finding pairs within distance threshold: %d bp", distance))

# Loop through signals and find pairs
for(i in 1:nrow(signals)) {
  for(j in (i+1):nrow(signals)) {
    if(j > nrow(signals)) break
    
    # Check if on same chromosome and within distance
    if(signals$chromosome[i] == signals$chromosome[j] &&
       abs(signals$position[i] - signals$position[j]) <= distance &&
       signals$gwas[i] != signals$gwas[j]) {  # Only compare different GWAS studies
      
      # Create unique pair ID
      pair_id <- paste0("pair_", i, "_", j)
      
      # Add to pairs dataframe
      pairs <- rbind(pairs, data.frame(
        pair_id=pair_id,
        gwas1=signals$gwas[i],
        gwas2=signals$gwas[j],
        chr=signals$chromosome[i],
        start=min(signals$position[i], signals$position[j]) - distance,
        end=max(signals$position[i], signals$position[j]) + distance
      ))
    }
  }
}

log_msg(sprintf("Found %d potential pairs", nrow(pairs)))

# Ensure start positions are not negative
pairs$start <- pmax(1, pairs$start)
log_msg("Adjusted start positions to ensure non-negative values")

# Write output
log_msg("Writing pairs output file")
write.csv(pairs, snakemake@output[["pairs"]], row.names=FALSE)
log_msg("Completed pair finding process")
