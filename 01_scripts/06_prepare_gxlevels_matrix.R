## Prepare input from eXpress to input into edgeR
## requires input of <SAMPLE>_L001_ilvd_trimmed.fq_non_rRNA.fq.sorted.bam_results.xprs
## 2018-11-27

# Clean space
# rm(list=ls())

# Details for output filename
date <- "2018-11-27"
ref.genome <- "16_libs_contig"
version <- "v.0.2"


# Set working directory
setwd("~/Documents/02_eRNA_taxo/eRNA_taxo/08_gx_levels/")

# Identify file names from folder in repo
files <- list.files(
  #path = "08_gx_levels", 
  pattern = "*_results.xprs")

# Read in required columns from files, put into a list format
current.name <- NULL;
expr.list <- list()

for(i in 1:length(files)){
  
  # obtain the name of the current file
  current.name <- files[i]
  print(current.name)
  
  # Reduce filename
  current.name <- gsub(x = current.name, pattern = "_L001_ilvd_trimmed.fq_non_rRNA.fq.sorted.bam_results.xprs",replacement = "")
  print(current.name)
  
  # Put the current file into a list, dropping the header row and keeping cols 2 (target ID) and 8 (eff_counts)
  expr.list[[i]] <- read.table(file = files[i], header = F, sep = "\t")[-1, c(2,8)]
  # Name the retained columns
  colnames(expr.list[[i]]) <- c(paste(current.name,".transcript", sep = ""),
                                paste(current.name,".eff.counts", sep = ""))
  
  # Sort the current list accession
  sort.column.name <- colnames(expr.list[[i]])[1]
  
  expr.list[[i]] <- 
    expr.list[[i]][with(expr.list[[i]], order(expr.list[[i]][sort.column.name])), ]
  
}

str(expr.list)

# this results in a list of all necessary info

# Just in case, you may save out the result, as this is a large object that can take a lot of time to produce
# save.image(file = "expr.list.RData")


# Build final matrix file
# Create variable with first column as the transcript ID
out.matrix <- expr.list[[1]][1]

# Build the matrix one sample at a time, collecting only the .eff.counts column
for(f in 1:length(expr.list)){
  out.matrix <- cbind(out.matrix, expr.list[[f]][2])
}

head(out.matrix)

# Rename the first column to reflect that it is the contig name
colnames(out.matrix)[1] <- "transcript.id"

# Save results
filename <- paste("out_matrix", "ref", ref.genome, version, date, sep = "_")
filename <- paste0(filename, ".csv")
filename 

write.csv(x = out.matrix, file = filename, row.names = F)
