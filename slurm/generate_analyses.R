library(readr)
library(dplyr)
library(tidyr)

signals <- read.csv("lead_snps.all_studies.csv")

signals$trait <- gsub(".*_.*_(.*)_.*", "\\1", signals$gwas)

window <- 250000
analyses <- data.frame()

for(i in 1:nrow(signals)){
  chr <- signals$chromosome[i]
  pos <- signals$position[i]
  
  matches <- which(signals$chromosome == chr &
                     abs(pos - signals$position) <= window)
  
  for(j in matches){
    if(signals$gwas[i] != signals$gwas[j]){
      tmp.df <- data.frame(gwas1 = signals$gwas[i], 
                           gwas2 = signals$gwas[j],
                           chr = chr,
                           pos.start = min(pos, signals$position[j]) - window,
                           pos.end = max(pos, signals$position[j]) + window)
      
      analyses <- rbind(analyses, tmp.df)
    }
  }
  
}

# adjust negative numbers
analyses$pos.start[analyses$pos.start < 0] <- 1

# remove duplicates
analyses$first <- apply(analyses[,1:2], 1, min)
analyses$second <- apply(analyses[,1:2], 1, max)
analyses2 <- analyses[,c(6,7,3:5)]
analyses2 <- unique(analyses2)
# extract traits
analyses$trait1 <- gsub(".*_.*_(.*)_.*", "\\1", analyses$gwas1)
analyses$trait2 <- gsub(".*_.*_(.*)_.*", "\\1", analyses$gwas2)

# add index
analyses2$analysis <- 1:nrow(analyses2)

#save
write_tsv(analyses2, "pairwise_analyses.txt")