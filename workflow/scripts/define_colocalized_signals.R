library(readr)
library(dplyr)
library(tidyr)

# Import data
coloc_results <- read_tsv(snakemake@input[["input_coloc_results"]], col_names=FALSE)
analysis_info <- read_csv(snakemake@input[["input_pairs"]])

colnames(coloc_results) <- c("analysis", "nvar", "h0", "h1", "h2", "h3", "h4", "min_p12")

# Join analysis ids
coloc_results <- left_join(coloc_results, analysis_info, by=c("analysis"="pair_index"))

# Add annotations
coloc_results$dense <- coloc_results$nvar >= as.numeric(snakemake@params[["dense_nvar_threshold"]])
coloc_results$dropout <- coloc_results$h0 >= as.numeric(snakemake@params[["signal_dropout_threshold"]])
                        coloc_results$h1 >= as.numeric(snakemake@params[["signal_dropout_threshold"]])
                        coloc_results$h2 >= as.numeric(snakemake@params[["signal_dropout_threshold"]])
coloc_results$coloc <- coloc_results$h4 >= as.numeric(snakemake@params[["coloc_h4_threshold"]])
coloc_results$robust <- coloc_results$min_p12 < as.numeric(snakemake@params[["sensitivity_p12_threshold"]])

# Check for AB/BA consistency 
coloc_consistent <- coloc_results %>%
    group_by(pair_id) %>%
    summarise(N = n(),
              dense = sum(dense),
              dropout = sum(dropout),
              coloc = sum(coloc),
              robust = sum(robust))

# Define colocalized signals
coloc_consistent <- coloc_consistent %>%
    mutate(colocalized = (N == as.numeric(snakemake@params[["min_n"]])) & 
                            (dense >= as.numeric(snakemake@params[["min_dense"]])) & 
                            (dropout <= as.numeric(snakemake@params[["max_dropout"]])) & 
                            (coloc >= as.numeric(snakemake@params[["min_coloc"]])) & 
                            (robust >= as.numeric(snakemake@params[["min_robust"]])))

# Create output
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

write_csv(output_df, snakemake@output[["output_coloc_signals"]])
write_csv(coloc_codes, snakemake@output[["output_coloc_codes"]])