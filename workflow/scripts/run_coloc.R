library(coloc)
library(readr)

# Read the extracted regions
combined_sumstats <- read_csv(snakemake@input[["input_joined_sumstats"]])
# Read input signals
signals <- read_csv(snakemake@input[["input_signals"]])
# Analysis data
pairwise_analyses <- read_csv(snakemake@input[["input_pairs"]])

# Filter to just columns of interest (ie drop allele columns)
wanted_cols <- c("chromosome", "position", "beta.gwas1", "se.gwas1", "beta.gwas2", "se.gwas2")
combined_sumstats <- combined_sumstats[, wanted_cols]

# QC
combined_sumstats[] <- sapply(combined_sumstats, as.numeric)
# If se for either GWAS is 0 or NA, remove variant
combined_sumstats <- combined_sumstats[combined_sumstats$se.gwas1 != 0 & combined_sumstats$se.gwas2 != 0,]
combined_sumstats <- combined_sumstats[!is.na(combined_sumstats$se.gwas1) & !is.na(combined_sumstats$se.gwas2),]

# Get study names and meta information
analysis_id <- snakemake@params[["pair_id"]]
analysis_meta <- as.data.frame(pairwise_analyses[pairwise_analyses$pair_index == analysis_id,])

if(nrow(analysis_meta) == 0) {
  stop("Pair ID not found in pairs file: ", analysis_id)
} else if(nrow(analysis_meta) > 1) {
  stop("Multiple entries found for Pair ID in pairs file: ", analysis_id)
}

gwas1_name <- analysis_meta$gwas1[1]
gwas2_name <- analysis_meta$gwas2[1]

gwas1_type <- as.character(signals$type[signals$gwas == gwas1_name][1])
gwas2_type <- as.character(signals$type[signals$gwas == gwas2_name][1])

gwas1_scalar <- as.numeric(signals$scalar[signals$gwas == gwas1_name][1])
gwas2_scalar <- as.numeric(signals$scalar[signals$gwas == gwas2_name][1])

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

# Run coloc
coloc_res <- coloc.abf(dataset1 = study1,
                       dataset2 = study2)


# Format results for output
output <- data.frame(
  nsnps = results$summary["nsnps"],
  PP.H0 = results$summary["PP.H0"],
  PP.H1 = results$summary["PP.H1"],
  PP.H2 = results$summary["PP.H2"],
  PP.H3 = results$summary["PP.H3"],
  PP.H4 = results$summary["PP.H4"]
)

# Write results
write_tsv(output, snakemake@output[["output_results"]])
