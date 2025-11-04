library(readr)
library(dplyr)

log_msg <- function(...) {
  msg <- paste(..., collapse = "")
  timestamped_msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  message(timestamped_msg)
  if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
    cat(timestamped_msg, file = snakemake@log[[1]], append = TRUE, sep = "\n")
  }
}

# Read in colocalized signals file
log_msg("Reading colocalized signals file")
coloc_signals <- read_csv(snakemake@input[["input_coloc_codes"]])
log_msg(sprintf("Loaded %d colocalized signals", nrow(coloc_signals)))

# Get info for this colocalized signal pair
log_msg(sprintf("Processing colocalized signal ID: %s", snakemake@params[["signal_code"]]))
coloc_signal <- coloc_signals[coloc_signals$coloc_id == snakemake@params[["signal_code"]], ]

if(nrow(coloc_signal) == 0) {
  stop("Colocalized signal ID not found in input file: ", snakemake@params[["signal_code"]])
}
log_msg("Found matching colocalized signal information")

# Load credible sets
log_msg(sprintf("Loading credible set for analysis1: %s", coloc_signal$analysis1))
cs1 <- read.delim(file.path(snakemake@input[["input_cs_dir"]], 
                          paste0(coloc_signal$analysis1, "_cs95.txt")), sep = "\t")
log_msg(sprintf("Loaded %d SNPs from analysis1 credible set", nrow(cs1)))

log_msg(sprintf("Loading credible set for analysis2: %s", coloc_signal$analysis2))
cs2 <- read.delim(file.path(snakemake@input[["input_cs_dir"]], 
                          paste0(coloc_signal$analysis2, "_cs95.txt")), sep = "\t")
log_msg(sprintf("Loaded %d SNPs from analysis2 credible set", nrow(cs2)))

# Join them
log_msg("Joining credible sets by chromosome and SNP")
combined_cs <- full_join(cs1, cs2, 
                          by=c("chr", "snp"), 
                          suffix=c(".analysis1", ".analysis2"))
log_msg(sprintf("Combined credible set has %d SNPs", nrow(combined_cs)))

# Save
log_msg("Writing joined credible sets output file")
write_csv(combined_cs, snakemake@output[["output_joined_cs"]])
log_msg("Completed joining credible sets")