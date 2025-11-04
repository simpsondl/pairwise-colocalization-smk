library(dplyr)
library(readr)
library(tidyr)

# Helper to log messages to both console and Snakemake log file (if available)
log_msg <- function(...) {
  msg <- paste(..., collapse = "")
  message(msg)
  if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
    cat(msg, file = snakemake@log[[1]], append = TRUE, sep = "\n")
  }
}

# Import data and settings
log_msg("Starting generate_analyses.R")
log_msg("Reading input signals from: ", snakemake@input[["input_signals"]])
signals <- read_csv(snakemake@input[["input_signals"]])
window <- snakemake@params[["distance"]]
log_msg("Distance window: ", window)
log_msg("Number of input signals: ", nrow(signals))

# Check if required columns are present in signals file
required_cols <- c("gwas", "chromosome", "position")
missing_cols <- setdiff(required_cols, colnames(signals))
if(length(missing_cols) > 0){
  log_msg("ERROR: Missing required columns in signals file: ", paste(missing_cols, collapse = ", "))
  stop(paste("Missing required columns in signals file:", paste(missing_cols, collapse = ", ")))
}

# Placeholders
analyses <- data.frame()
exact_matches <- data.frame()

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
      # Exclude same SNP comparisons
      if(pos != signals$position[j]){
        # Define region of interest around the two nearby signals  
        tmp.df <- data.frame(gwas1 = signals$gwas[i], 
                             gwas2 = signals$gwas[j],
                             chr = chr,
                             gwas1_signal_pos = pos,
                             gwas2_signal_pos = signals$position[j],
                             pos.start = pos - window,
                             pos.end = pos + window)
        analyses <- rbind(analyses, tmp.df)
      }  else {
        # If same variant is identified, no need to colocalize
        # just record the signal position
        tmp.df <- data.frame(gwas1 = signals$gwas[i], 
                             gwas2 = signals$gwas[j],
                             chr = chr,
                             gwas_signal_pos = pos)
        exact_matches <- rbind(exact_matches, tmp.df)
      }    
    }
  } 
}

# adjust negative numbers
if(nrow(analyses) > 0){
  analyses$pos.start[analyses$pos.start < 0] <- 1
}

# add id to represent an analysis up to order
analyses <- analyses %>%
  rowwise() %>%
  mutate(pair_id = paste(paste(sort(c(gwas1_signal_pos, gwas2_signal_pos)), collapse = "_"),
                         paste(sort(c(gwas1, gwas2)), collapse = "_"), sep = "_"), .before = 1) %>%
  ungroup() %>%
  arrange(pair_id)

# add a unique index for each analysis
if(nrow(analyses) > 0){
  analyses <- analyses %>% 
    mutate(pair_index = paste0("pw_",1:nrow(analyses)), .before = 1) 
} else {
  analyses <- analyses %>% mutate(pair_index = character(0))
}

# save
log_msg("Number of analyses generated: ", nrow(analyses))
log_msg("Number of exact matches: ", nrow(exact_matches))
log_msg("Writing analyses to: ", snakemake@output[["output_pairs"]])
write_csv(analyses, snakemake@output[["output_pairs"]])
log_msg("Writing exact matches to: ", snakemake@output[["output_exact_matches"]])
write_csv(exact_matches, snakemake@output[["output_exact_matches"]])
log_msg("generate_analyses.R finished successfully")