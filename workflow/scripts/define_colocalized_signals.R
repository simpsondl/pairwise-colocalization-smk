library(readr)
library(dplyr)
library(tidyr)

log_msg <- function(...) {
  msg <- paste(..., collapse = "")
  timestamped_msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  message(timestamped_msg)
  if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
    cat(timestamped_msg, file = snakemake@log[[1]], append = TRUE, sep = "\n")
  }
}

# Import data
log_msg("Reading colocalization results file")
coloc_results <- read_tsv(snakemake@input[["input_coloc_results"]], col_names=FALSE)
log_msg(sprintf("Loaded %d colocalization results", nrow(coloc_results)))

log_msg("Reading analysis pairs file")
analysis_info <- read_csv(snakemake@input[["input_pairs"]])
log_msg(sprintf("Loaded %d analysis pairs", nrow(analysis_info)))

colnames(coloc_results) <- c("analysis", "nvar", "h0", "h1", "h2", "h3", "h4", "min_p12")

# Join analysis ids
log_msg("Joining colocalization results with analysis information")
coloc_results <- left_join(coloc_results, analysis_info, by=c("analysis"="pair_index"))
log_msg(sprintf("Joined data has %d rows", nrow(coloc_results)))

# Add annotations
log_msg("Adding quality control annotations")
coloc_results$dense <- coloc_results$nvar >= as.numeric(snakemake@params[["dense_nvar_threshold"]])
coloc_results$dropout <- coloc_results$h0 >= as.numeric(snakemake@params[["signal_dropout_threshold"]])
                        coloc_results$h1 >= as.numeric(snakemake@params[["signal_dropout_threshold"]])
                        coloc_results$h2 >= as.numeric(snakemake@params[["signal_dropout_threshold"]])
coloc_results$coloc <- coloc_results$h4 >= as.numeric(snakemake@params[["coloc_h4_threshold"]])
coloc_results$robust <- coloc_results$min_p12 < as.numeric(snakemake@params[["sensitivity_p12_threshold"]])

# Check for AB/BA consistency 
log_msg("Checking for AB/BA consistency in colocalization results")
coloc_consistent <- coloc_results %>%
    group_by(pair_id) %>%
    summarise(N = n(),
              dense = sum(dense),
              dropout = sum(dropout),
              coloc = sum(coloc),
              robust = sum(robust))

log_msg(sprintf("Found %d unique pairs for consistency check", nrow(coloc_consistent)))

# Define colocalized signals
log_msg("Defining colocalized signals based on consistency criteria")
coloc_consistent <- coloc_consistent %>%
    mutate(colocalized = (N == as.numeric(snakemake@params[["min_n"]])) & 
                            (dense >= as.numeric(snakemake@params[["min_dense"]])) & 
                            (dropout <= as.numeric(snakemake@params[["max_dropout"]])) & 
                            (coloc >= as.numeric(snakemake@params[["min_coloc"]])) & 
                            (robust >= as.numeric(snakemake@params[["min_robust"]])))

log_msg(sprintf("Identified %d colocalized signal pairs", sum(coloc_consistent$colocalized)))

# Create output
log_msg("Creating output dataframes")
output_df <- analysis_info %>%
    filter(pair_id %in% coloc_consistent$pair_id[coloc_consistent$colocalized]) %>%
    distinct(pair_id, .keep_all=TRUE) %>%
    select(gwas1, gwas2, chr, gwas1_signal_pos, gwas2_signal_pos)

coloc_codes <- coloc_results %>%
    filter(pair_id %in% coloc_consistent$pair_id[coloc_consistent$colocalized]) %>%
    select(pair_id, analysis) %>%
    group_by(pair_id) %>%
    summarise(analysis1 = analysis[1],
              analysis2 = analysis[2]) %>%
    mutate(coloc_id = paste0("coloc_",1:nrow(.)))

log_msg("Writing colocalized signals output file")
write_csv(output_df, snakemake@output[["output_coloc_signals"]])
log_msg("Writing colocalization codes output file")
write_csv(coloc_codes, snakemake@output[["output_coloc_codes"]])
log_msg("Completed defining colocalized signals")