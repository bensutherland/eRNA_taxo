## Low expression filter and normalization of count data from eXpress output
# Will produce "08_gx_levels/expr_contigs.csv"

# Clear space
#rm(list=ls())

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
rownames(my.counts) <- my.counts[,1]
my.counts.round <- round(my.counts[,-1])
str(my.counts)
str(my.counts.round)

# create DGElist
my.counts <- DGEList(counts = my.counts.round)

#### 2. Filter Data ####
# Find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.reads.mapping.per.transcript <- 5
cpm.filt <- min.reads.mapping.per.transcript / min(my.counts$samples$lib.size) * 1000000
cpm.filt # min cpm filt

min.ind <- 5 # choose the minimum number of individuals that need to pass the threshold

# identify tags passing filter
keep <- rowSums(cpm(my.counts)>cpm.filt) >= min.ind # Find which transcripts pass the filter
table(keep) # gives number passing, number failing
##16_lib: total retained (TRUE) = 32,866 transcripts of 7,971,030
##bbmap_only: total retained (TRUE) =    
##cdhit95_bbmap: total retained (TRUE) Logan = 20488 # Xavier 20363

# subset DGEList
my.counts <- my.counts[keep, , keep.lib.sizes=FALSE] #keep.lib.sizes = T retains original lib sizes, otherwise recomputes w remaining tags
dim(my.counts)

#### 3. Normalization ####
# Use TMM normalization, as it takes into account highly expressed genes that may take up sequencing rxn and make other genes look down-reg.
my.counts <- calcNormFactors(my.counts, method = c("TMM"))
my.counts$samples

# Plot norm.factors by library size
plot(my.counts$samples$norm.factors ~ my.counts$samples$lib.size)


#### 4. Save results and export background list for annotation ####
# save normalized image
save.image(file = "08_gx_levels/normalized.RData")

# Export background list for annotation
expr.contigs <- dimnames(my.counts)[[1]]
write.table(x = expr.contigs, sep = ",", file = "08_gx_levels/expr_contigs.csv", row.names = F, quote = F, col.names = F)

