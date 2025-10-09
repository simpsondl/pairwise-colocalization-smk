library(dplyr)

# Read pairs file
pairs <- read.csv(snakemake@input[["pairs"]])

# Get info for this pair
pair <- pairs[pairs$pair_index == snakemake@params[["pair_id"]], ]

if(nrow(pair) == 0) {
  stop("Pair ID not found in pairs file: ", snakemake@params[["pair_id"]])
}

# Function to extract region from a GWAS file
extract_region <- function(gwas_file, chr, start, end) {
  # Read the full file 
  data <- read.table(gwas_file, header=FALSE)
  colnames(data) <- c("chromosome",	"position",	"allele1",	"allele2",	"beta",	"se",	"pval")
  
  # Extract region
  region_data <- data[data$chromosome == chr & 
                     data$position >= start & 
                     data$position <= end, ]
  
  return(region_data)
}

# Extract regions for both GWAS
gwas1 <- extract_region(
  file.path(snakemake@input[["sumstats_dir"]], paste0(pair$gwas1, ".txt")),
  pair$chr, pair$pos.start, pair$pos.end
)

gwas2<- extract_region(
  file.path(snakemake@input[["sumstats_dir"]], paste0(pair$gwas2, ".txt")),
  pair$chr, pair$pos.start, pair$pos.end
)

both_gwas <- inner_join(gwas1[,1:6], gwas2[,1:6], 
                        by=c("chromosome", "position", "allele1", "allele2"), 
                        suffix=c(".gwas1", ".gwas2"))

write.csv(both_gwas, snakemake@output[["output_joined_sumstats"]], row.names=FALSE, quote=FALSE)