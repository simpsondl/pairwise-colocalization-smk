library(readr)
library(coloc)

# Get filepath to sumstats file from argument
args <- commandArgs(trailingOnly = TRUE)
sumstats_fpath <- args[[1]]
print(sumstats_fpath)

# Parse analysis id
analysis_id <- as.numeric(sub("^[^_]*_([^_]+)_.*", "\\1", 
                              gsub(".*/", "", sumstats_fpath)))

# Import data
study_meta <- read.delim("gwas_study_meta_information.txt")
pairwise_coloc <- read.delim("pairwise_coloc.txt", header = FALSE)
combined_sumstats <- read_table(sumstats_fpath, col_names = FALSE)

analysis_meta <- as.data.frame(pairwise_coloc[pairwise_coloc$V6 == analysis_id,])

# QC
combined_sumstats[] <- sapply(combined_sumstats, as.numeric)
combined_sumstats <- combined_sumstats[combined_sumstats$X3 != 0 & combined_sumstats$X5 != 0,]
combined_sumstats <- combined_sumstats[!is.na(combined_sumstats$X3) & !is.na(combined_sumstats$X5),]

# Set up study data structures
study1 <- list(beta = combined_sumstats$X2,
               varbeta = (combined_sumstats$X3)^2,
               snp = as.character(combined_sumstats$X1),
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

analysis_results <- data.frame(analysis = analysis_id,
			       nvar = nrow(combined_sumstats),
                               h0 = coloc_res$summary[2],
                               h1 = coloc_res$summary[3],
                               h2 = coloc_res$summary[4],
                               h3 = coloc_res$summary[5],
                               h4 = coloc_res$summary[6])

if(analysis_results$h4 > 0.5){
  # Make a credible set
  o <- order(coloc_res$results$SNP.PP.H4, decreasing = TRUE)
  cs <- cumsum(coloc_res$results$SNP.PP.H4[o])
  w <- which(cs > 0.95)[1]
  
  results_output <- coloc_res$results[o,][1:w,][,c(1,11)]
  results_output$chr <- unique(analysis_meta$V3) 
  results_output <- results_output[,c(3,1,2)]
}

write_tsv(analysis_results, 
          paste0("coloc_outputs/analysis_", analysis_id, "_results.txt"))

if(exists("results_output")){
  write_tsv(results_output, 
            paste0("coloc_outputs/analysis_", analysis_id, "_credibleset.txt"))  
}
