library(coloc)
library(readr)

# Logging helper
log_msg <- function(...) {
  msg <- paste(..., collapse = "")
  timestamped_msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  message(timestamped_msg)
  if (!is.null(snakemake@log) && length(snakemake@log) >= 1) {
    cat(timestamped_msg, file = snakemake@log[[1]], append = TRUE, sep = "\n")
  }
}

log_msg("Starting run_coloc.R for pair_id: ", snakemake@params[["pair_id"]])

# Read the extracted regions
log_msg("Reading combined sumstats from: ", snakemake@input[["input_joined_sumstats"]])
combined_sumstats <- read_csv(snakemake@input[["input_joined_sumstats"]])
# Read input signals
log_msg("Reading input signals from: ", snakemake@input[["input_signals"]])
signals <- read_csv(snakemake@input[["input_signals"]])
# Analysis data
log_msg("Reading pairwise analyses from: ", snakemake@input[["input_pairs"]])
pairwise_analyses <- read_csv(snakemake@input[["input_pairs"]])

log_msg("Initial combined sumstats has ", nrow(combined_sumstats), " rows")

# Filter to just columns of interest (ie drop allele columns)
wanted_cols <- c("chromosome", "position", "beta.gwas1", "se.gwas1", "beta.gwas2", "se.gwas2")
combined_sumstats <- combined_sumstats[, wanted_cols]

# QC
combined_sumstats[] <- sapply(combined_sumstats, as.numeric)
# If se for either GWAS is 0 or NA, remove variant
combined_sumstats <- combined_sumstats[combined_sumstats$se.gwas1 != 0 & combined_sumstats$se.gwas2 != 0,]
combined_sumstats <- combined_sumstats[!is.na(combined_sumstats$se.gwas1) & !is.na(combined_sumstats$se.gwas2),]

log_msg("After QC, combined sumstats has ", nrow(combined_sumstats), " rows")

# Get study names and meta information
analysis_id <- snakemake@params[["pair_id"]]
analysis_meta <- as.data.frame(pairwise_analyses[pairwise_analyses$pair_index == analysis_id,])
print(colnames(analysis_meta))

if(nrow(analysis_meta) == 0) {
  log_msg("ERROR: Pair ID not found in pairs file: ", analysis_id)
  stop("Pair ID not found in pairs file: ", analysis_id)
} else if(nrow(analysis_meta) > 1) {
  log_msg("ERROR: Multiple entries found for Pair ID in pairs file: ", analysis_id)
  stop("Multiple entries found for Pair ID in pairs file: ", analysis_id)
}

gwas1_name <- analysis_meta$gwas1[1]
gwas2_name <- analysis_meta$gwas2[1]

gwas1_type <- as.character(signals$type[signals$gwas == gwas1_name][1])
gwas2_type <- as.character(signals$type[signals$gwas == gwas2_name][1])

gwas1_scalar <- as.numeric(signals$scalar[signals$gwas == gwas1_name][1])
gwas2_scalar <- as.numeric(signals$scalar[signals$gwas == gwas2_name][1])

log_msg("GWAS1: ", gwas1_name, " (type: ", gwas1_type, ", scalar: ", gwas1_scalar, ")")
log_msg("GWAS2: ", gwas2_name, " (type: ", gwas2_type, ", scalar: ", gwas2_scalar, ")")

# Set up study data structures
study1 <- list(beta = combined_sumstats$beta.gwas1,
               varbeta = (combined_sumstats$se.gwas1)^2,
               snp = paste(combined_sumstats$chromosome, combined_sumstats$position, sep = ":"),
               type = gwas1_type,
               scalar = gwas1_scalar)

ifelse(study1$type == "cc", 
       names(study1)[names(study1) == "scalar"] <- "s", 
       names(study1)[names(study1) == "scalar"] <- "sdY")

study2 <- list(beta = combined_sumstats$beta.gwas2,
               varbeta = (combined_sumstats$se.gwas2)^2,
               snp = paste(combined_sumstats$chromosome, combined_sumstats$position, sep = ":"),
               type = gwas2_type,
               scalar = gwas2_scalar)

ifelse(study2$type == "cc", 
       names(study2)[names(study2) == "scalar"] <- "s", 
       names(study2)[names(study2) == "scalar"] <- "sdY")

log_msg("Study data structures set up with ", length(study1$snp), " SNPs")

# Run coloc
log_msg("Running coloc.abf...")
coloc_res <- coloc.abf(dataset1 = study1,
                       dataset2 = study2)

log_msg("Coloc results: H0=", coloc_res$summary[2], ", H1=", coloc_res$summary[3], 
         ", H2=", coloc_res$summary[4], ", H3=", coloc_res$summary[5], ", H4=", coloc_res$summary[6])

# Perform sensitivity analysis
log_msg("Running sensitivity analysis...")
coloc_sens <- sensitivity(coloc_res, rule="H4 > 0.5", 
                          doplot = FALSE, plot.manhattans = FALSE)

log_msg("Sensitivity analysis completed")

# Format results for output
analysis_results <- data.frame(analysis = analysis_id,
			                         nvar = nrow(combined_sumstats),
                               h0 = coloc_res$summary[2],
                               h1 = coloc_res$summary[3],
                               h2 = coloc_res$summary[4],
                               h3 = coloc_res$summary[5],
                               h4 = coloc_res$summary[6],
                               min_p12 = min(coloc_sens$p12[coloc_sens$pass]))

# Check if we should make a credible set
if(analysis_results$h4 > 0.5){
  log_msg("H4 > 0.5, creating credible set...")
  # Make a credible set
  o <- order(coloc_res$results$SNP.PP.H4, decreasing = TRUE)
  cs <- cumsum(coloc_res$results$SNP.PP.H4[o])
  w <- which(cs > 0.95)[1]
  
  cs_output <- coloc_res$results[o,][1:w,][,c(1,11)]
  cs_output$chr <- unique(analysis_meta$chr) 
  cs_output <- cs_output[,c(3,1,2)]
  log_msg("Credible set created with ", nrow(cs_output), " SNPs")
} else {
  cs_output <- data.frame()
  log_msg("H4 <= 0.5, no credible set created")
}

# Write results
log_msg("Writing results to: ", snakemake@output[["output_results"]])
write_tsv(analysis_results, snakemake@output[["output_results"]])
log_msg("Writing credible set to: ", snakemake@output[["output_cs"]])
write_tsv(cs_output, snakemake@output[["output_cs"]])

log_msg("run_coloc.R finished successfully for pair_id: ", analysis_id)

