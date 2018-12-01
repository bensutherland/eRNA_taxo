## Low expression filter and normalization of count data from eXpress output
# Will produce "08_gx_levels/expr_contigs.csv"

# Clear space
# rm(list=ls())

## Install packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
require("edgeR")
#install.packages("readr")
require("readr")

#biocLite("locfit")
require("locfit")

# User guides are here:
# browseVignettes("edgeR")
# edgeRUsersGuide()

# Set working directory
setwd("~/Documents/02_eRNA_taxo/eRNA_taxo")

#### 1. Import Data ####
# Import interpretation file
interp <- as.data.frame(read.csv("00_archive/eRNA_interp_2017-09-05.csv"))
head(interp)
names(interp)

# Import counts file (takes it out of 08_gx_levels)
out.matrix.filename <- list.files(path = "08_gx_levels/", pattern = "out_matrix_*")
out.matrix.filename <- paste0("08_gx_levels/", out.matrix.filename)

my.counts <- read.csv(file = out.matrix.filename)

# Import data
dim(my.counts) 

## Backup stuff
# my.counts.bk <- my.counts # just in case
# my.counts <- my.counts.bk # CAUTION go backwarks

# Set up DGEList
rownames(my.counts) <- my.counts[,1] # contig names assigned as rownames
my.counts.round <- round(my.counts[,-1]) # round expr values to nearest whole number
str(my.counts)
str(my.counts.round) # to confirm

# Create DGElist
my.counts <- DGEList(counts = my.counts.round)

#### 2. Filter Data Based on Low Expression ####
# Find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.reads.mapping.per.transcript <- 5
cpm.filt <- min.reads.mapping.per.transcript / min(my.counts$samples$lib.size) * 1000000
cpm.filt # min cpm filt

# Minimum number of individuals needed to call the transcript expressed
min.ind <- 5

# Identify tags passing filter
keep <- rowSums(cpm(my.counts)>cpm.filt) >= min.ind # Find which transcripts pass the filter
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
