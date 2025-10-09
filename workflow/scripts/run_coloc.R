library(coloc)
library(readr)

# Read the extracted regions
combined_sumstats <- read_csv(snakemake@input[["input_joined_sumstats"]])
# Read input signals
signals <- read_csv(snakemake@input[["input_signals"]])



# QC
combined_sumstats[] <- sapply(combined_sumstats, as.numeric)
combined_sumstats <- combined_sumstats[combined_sumstats$se.gwas1 != 0 & combined_sumstats$se.gwas2 != 0,]
combined_sumstats <- combined_sumstats[!is.na(combined_sumstats$se.gwas1) & !is.na(combined_sumstats$se.gwas2),]

# Set up study data structures
study1 <- list(beta = combined_sumstats$beta.gwas1,
               varbeta = (combined_sumstats$se.gwas1)^2,
               snp = paste(combined_sumstats$chromosome, combined_sumstats$position, sep = ":"),
               type = study_meta$Type[study_meta$Data_File == analysis_meta$V1[1]])

if(study1$type == "cc"){
  study1 <- c(study1, 
                 s = study_meta$s[study_meta$Data_File == analysis_meta$V1[1]])
} else {
  study1 <- c(study1,
                 sdY = 1)
}

study2 <- list(beta = combined_sumstats$X4,
               varbeta = (combined_sumstats$X5)^2,
               snp = as.character(combined_sumstats$X1),
               type = study_meta$Type[study_meta$Data_File == analysis_meta$V2[1]])

if(study2$type == "cc"){
  study2 <- c(study2, 
                 s = study_meta$s[study_meta$Data_File == analysis_meta$V2[1]])
} else {
  study2 <- c(study2,
                 sdY = 1)
}

print(analysis_id)
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
write.table(output, snakemake@output[["results"]], 
            row.names=FALSE, quote=FALSE, sep="\t")
