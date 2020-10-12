## Low expression filter and normalization of count data from eXpress output
# Will produce "08_gx_levels/expr_contigs.csv"

# Clear space
# rm(list=ls())

## Install packages
#install.packages("BiocManager")
packageVersion(pkg = "BiocManager") # v.1.30.10
#BiocManager::install("edgeR")
require("edgeR")
packageVersion(pkg = "edgeR") # v.3.28.1

#install.packages("readr")
require("readr")
packageVersion(pkg = "readr") # v.1.4.0

#BiocManager::install("locfit")
require("locfit")
packageVersion(pkg = "locfit") #v.1.5.9.4

# Set working directory
setwd("/hdd/20_other_research/metatranscriptomics/eRNA_taxo")


#### 1. Import Data ####
# Import interpretation file
interp <- as.data.frame(read.csv("00_archive/eRNA_interp_2017-09-05.csv"))
head(interp)
dim(interp)
names(interp)

# Import read count file from 08_gx_levels (requires only one out_matrix in target folder)
# Identify file
out.matrix.filename <- list.files(path = "08_gx_levels/", pattern = "out_matrix_*")
out.matrix.filename <- paste0("08_gx_levels/", out.matrix.filename)
# Import file
my.counts <- read.csv(file = out.matrix.filename)

# View data
dim(my.counts)
my.counts[1:5,1:5]
max(my.counts[,c(2:ncol(my.counts))]) # max value is 49,482.88 (therefore is linear vals at this point)
sum(my.counts$S13_S3.eff.counts) # 1,296,760

## Backup stuff
# my.counts.bk <- my.counts # just in case
# my.counts <- my.counts.bk # CAUTION go backwarks

# Set up DGEList using linear, whole numbers
rownames(my.counts) <- my.counts[,1] # assign contig names as rownames
my.counts.round <- round(my.counts[,-1]) # round the linear effective count values to nearest whole number
str(my.counts)
str(my.counts.round) # to confirm that the function worked as expected

# Create DGElist
my.counts <- DGEList(counts = my.counts.round)


#### 2. Filter Data Based on Low Expression ####
# Find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.reads.mapping.per.transcript <- 5
cpm.filt <- min.reads.mapping.per.transcript / min(my.counts$samples$lib.size) * 1000000
cpm.filt
min.reads.mapping.per.transcript / max(my.counts$samples$lib.size) * 1000000 # just for comparison
# so we will consider a gene expressed if it has 3.86 alignments in a library size of 1,295,760
cpm.filt # min cpm filt

# Minimum number of individuals needed to call the transcript expressed
min.ind <- 5

# Identify tags passing filter
keep <- rowSums(cpm(my.counts) > cpm.filt) >= min.ind # Find which transcripts pass the filter
table(keep) # TRUE = expressed

# Subset DGEList based on these low expressed tags
my.counts <- my.counts[keep, , keep.lib.sizes=FALSE] # keep.lib.sizes = T retains original lib sizes, otherwise recomputes w remaining tags
dim(my.counts)


#### 3. Normalization ####
# Use TMM normalization, as it takes into account highly expressed genes that may take up sequencing rxn and make other genes look down-reg.
my.counts <- calcNormFactors(my.counts, method = c("TMM"))
my.counts$samples
plot(my.counts$samples$norm.factors ~ my.counts$samples$lib.size) # Plot norm.factors by library size


#### 5. Prepare Output ####
# Generate filtered, CPM matrix, linear or log2
normalized.output.linear <- cpm(my.counts, normalized.lib.sizes = T, log = F)
write.csv(normalized.output.linear, file = "08_gx_levels/normalized_output_linear.csv", quote = F)

normalized.output.log2 <- cpm(my.counts, normalized.lib.sizes = T, log = T, prior.count = 1)
write.csv(normalized.output.log2, file = "08_gx_levels/normalized_output_log2.csv", quote = F)

# Export background list for annotation
expr.contigs <- dimnames(my.counts)[[1]]
write.table(x = expr.contigs, sep = ",", file = "08_gx_levels/expr_contigs.csv", row.names = F, quote = F, col.names = F)

# Save overall normalized image
save.image(file = "08_gx_levels/normalized.RData")
