# Read input signals
signals <- read.csv(snakemake@input[["signals"]])
pairs <- data.frame()

# Find pairs within distance threshold
distance <- snakemake@params[["distance"]]

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

# Ensure start positions are not negative
pairs$start <- pmax(1, pairs$start)

# Write output
write.csv(pairs, snakemake@output[["pairs"]], row.names=FALSE)
