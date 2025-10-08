library(dplyr)
library(readr)
library(tidyr)

signals <- read_csv(snakemake@input[["input_signals"]])

write_csv(signals, snakemake@output[["output_pairs"]])