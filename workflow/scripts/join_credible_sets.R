library(readr)
library(dplyr)

# Read in colocalized signals file
coloc_signals <- read_csv(snakemake@input[["input_coloc_codes"]])

# Get info for this colocalized signal pair
coloc_signal <- coloc_signals[coloc_signals$coloc_id == snakemake@params[["signal_code"]], ]

if(nrow(coloc_signal) == 0) {
  stop("Colocalized signal ID not found in input file: ", snakemake@params[["signal_code"]])
}

# Load credible sets
cs1 <- read.delim(file.path(snakemake@input[["input_cs_dir"]], 
                          paste0(coloc_signal$analysis1, "_cs95.txt")), sep = "\t")

cs2 <- read.delim(file.path(snakemake@input[["input_cs_dir"]], 
                          paste0(coloc_signal$analysis2, "_cs95.txt")), sep = "\t")

# Join them
combined_cs <- full_join(cs1, cs2, 
                          by=c("chr", "snp"), 
                          suffix=c(".analysis1", ".analysis2"))

# Save
write_csv(combined_cs, snakemake@output[["output_joined_cs"]])