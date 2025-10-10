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
print(colnames(analysis_meta))

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
analysis_results <- data.frame(analysis = analysis_id,
			                         nvar = nrow(combined_sumstats),
                               h0 = coloc_res$summary[2],
                               h1 = coloc_res$summary[3],
                               h2 = coloc_res$summary[4],
                               h3 = coloc_res$summary[5],
                               h4 = coloc_res$summary[6])

# Check if we should make a credible set
if(analysis_results$h4 > 0.5){
  # Make a credible set
  o <- order(coloc_res$results$SNP.PP.H4, decreasing = TRUE)
  cs <- cumsum(coloc_res$results$SNP.PP.H4[o])
  w <- which(cs > 0.95)[1]
  
  cs_output <- coloc_res$results[o,][1:w,][,c(1,11)]
  cs_output$chr <- unique(analysis_meta$chr) 
  cs_output <- cs_output[,c(3,1,2)]
} else {
  cs_output <- data.frame()
}

# Write results
write_tsv(analysis_results, snakemake@output[["output_results"]])
write_tsv(cs_output, snakemake@output[["output_cs"]])

