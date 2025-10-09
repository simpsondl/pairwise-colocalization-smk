library(coloc)

# Read the extracted regions
gwas1 <- read.table(snakemake@input[["gwas1"]], header=TRUE)
gwas2 <- read.table(snakemake@input[["gwas2"]], header=TRUE)

# Set up dataset 1
dataset1 <- list(
  beta = gwas1$beta,
  varbeta = gwas1$se^2,
  type = "quant",
  snp = gwas1$rsid,
  position = gwas1$position
)

# Set up dataset 2
dataset2 <- list(
  beta = gwas2$beta,
  varbeta = gwas2$se^2,
  type = "quant",
  snp = gwas2$rsid,
  position = gwas2$position
)

# Run coloc
results <- coloc.abf(dataset1, dataset2)

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
