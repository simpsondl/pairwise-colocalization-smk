# Read pairs file
pairs <- read.csv(snakemake@input[["pairs"]])

# Get info for this pair
pair <- pairs[pairs$pair_index == snakemake@params[["pair_id"]], ]

if(nrow(pair) == 0) {
  stop("Pair ID not found in pairs file: ", snakemake@params[["pair_id"]])
}

# Function to extract region from a GWAS file
extract_region <- function(gwas_file, chr, start, end, out_file) {
  # Read the full file first to ensure proper column handling
  data <- read.table(gwas_file, header=FALSE)
  colnames(data) <- c("chromosome",	"position",	"allele1",	"allele2",	"beta",	"se",	"pval")
  
  # Extract region
  region_data <- data[data$chromosome == chr & 
                     data$position >= start & 
                     data$position <= end, ]
  
  # Write extracted region
  write.table(region_data, out_file, 
              row.names=FALSE, quote=FALSE, sep="\t")
}

# Extract regions for both GWAS
extract_region(
  file.path(snakemake@input[["sumstats_dir"]], paste0(pair$gwas1, ".txt")),
  pair$chr, pair$start, pair$end,
  snakemake@output[["gwas1"]]
)

extract_region(
  file.path(snakemake@input[["sumstats_dir"]], paste0(pair$gwas2, ".txt")),
  pair$chr, pair$start, pair$end,
  snakemake@output[["gwas2"]]
)
