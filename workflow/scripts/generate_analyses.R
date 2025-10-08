library(dplyr)
library(readr)
library(tidyr)

# Import data and settings
signals <- read_csv(snakemake@input[["input_signals"]])
window <- snakemake@config[["DISTANCE_LIMIT"]]

# Check if required columns are present in signals file
required_cols <- c("gwas", "chromosome", "position")
missing_cols <- setdiff(required_cols, colnames(signals))
if(length(missing_cols) > 0){
  stop(paste("Missing required columns in signals file:", paste(missing_cols, collapse = ", ")))
}

# Placeholder
analyses <- data.frame()

# Loop over each signal
for(i in 1:nrow(signals)){
  chr <- signals$chromosome[i]
  pos <- signals$position[i]
  
  matches <- which(signals$chromosome == chr &
                     abs(pos - signals$position) <= window)
  
  # Loop over all signals close to the current signal
  for(j in matches){
    # Exclude same GWAS comparisons
    if(signals$gwas[i] != signals$gwas[j]){
      # Define region of interest around the two nearby signals  
      tmp.df <- data.frame(gwas1 = signals$gwas[i], 
                           gwas2 = signals$gwas[j],
                           chr = chr,
                           pos.start = pos - window,
                           pos.end = pos + window)
      analyses <- rbind(analyses, tmp.df)
    }
  } 
}

# adjust negative numbers
analyses$pos.start[analyses$pos.start < 0] <- 1

# add an index
analyses <- analyses %>% mutate(pair_index = paste0("pw_",1:nrow(analyses)), .before = 1) 

# save
write_csv(analyses, snakemake@output[["output_pairs"]])