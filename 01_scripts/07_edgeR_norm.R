## Normalization of count data from eXpress output (in table)

#rm(list=ls())

## Install Packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
require("edgeR")

biocLite("locfit")
require("locfit")

# User guides are here:
# browseVignettes("edgeR")
# edgeRUsersGuide()

setwd("~/Documents/02_eRNA_oyster/eRNA_taxo")

#### 1. Import Data ####
# Import interpretation file
interp <- as.data.frame(read.csv("00_archive/eRNA_interp_2017-09-05.csv"))
head(interp)
names(interp)

# Import counts file
my.counts <- read.csv(file = "08_gx_levels/out.matrix.csv")

# Import data
dim(my.counts) #total unique tags = 4,528,162

## Backup stuff
# my.counts.bk <- my.counts # just in case
# my.counts <- my.counts.bk # go backwarks

# Set up DGEList
rownames(my.counts) <- my.counts[,1]
head(my.counts)

my.counts.round <- round(my.counts[,-1])
str(my.counts.round)

# create DGElist
my.counts <- DGEList(counts = my.counts.round)
# my.counts <- DGEList(counts = my.counts.round, group =) # can also create a group function for use in an exact test

#### 2. Filter Data ####
# Find an optimal cpm filt (edgeRuserguide suggests 5-10 reads mapping to transcript)
min.reads.mapping.per.transcript <- 5
cpm.filt <- min.reads.mapping.per.transcript / min(my.counts$samples$lib.size) * 1000000
cpm.filt # min cpm filt

min.ind <- 5 # choose the minimum number of individuals that need to pass the threshold

# identify tags passing filter
keep <- rowSums(cpm(my.counts)>cpm.filt) >= min.ind # Find which transcripts pass the filter
table(keep) # gives number passing, number failing

# subset DGEList
my.counts <- my.counts[keep, , keep.lib.sizes=FALSE] #keep.lib.sizes = T retains original lib sizes, otherwise recomputes w remaining tags
dim(my.counts)

#### 3. Normalization ####
# Use TMM normalization, as it takes into account highly expressed genes that may take up sequencing rxn and make other genes look down-reg.
my.counts <- calcNormFactors(my.counts, method = c("TMM"))
my.counts$samples
plot(my.counts$samples$norm.factors ~ my.counts$samples$lib.size)


#### 3. Create design matrix ####
interp$file.name!="S13"
#sampleType <- interp$Range
binary.pCO2 <- interp$Range[interp$file.name!="S13"] # currently missing one
binary.season <- interp$season[interp$file.name!="S13"] # currently missing one

designMat <- model.matrix(~binary.pCO2)
designMat <- model.matrix(~binary.pCO2 * binary.season)
designMat <- model.matrix(~binary.season)

# Estimate dispersions (measure inter-library variation per tag)
#my.counts <- estimateDisp(my.counts) # note that this can use a design matrix when provided 
my.counts <- estimateDisp(my.counts, design=designMat) # note that this can use a design matrix when provided 

my.counts <- estimateGLMCommonDisp(my.counts, design=designMat)
my.counts <- estimateGLMTrendedDisp(my.counts, design=designMat)
my.counts <- estimateGLMTagwiseDisp(my.counts, design=designMat)

summary(my.counts$prior.df) # est. overall var. across genome for dataset
sqrt(my.counts$common.disp) #coeff of var, for biol. var

plotBCV(my.counts)


# #### 4. Prepare Output ####
# # generate CPM matrix
# normalized.output <- cpm(my.counts, normalized.lib.sizes = TRUE, log= F)
# 
# # Compare the raw counts to the normalized cpm values (not log)
# my.counts$counts[1:5, 1:5] # not normalized, raw counts
# normalized.output[1:5, 1:5] # normalized lib size calculated cpm values
# 
# # # output as normalized linear
# # write.csv(normalized.output, file = "03_normalized_data/normalized_output_matrix.csv")
# # 
# # # # output as normalized log2 (in progress)
# # normalized.output.log2 <- cpm(my.counts, normalized.lib.sizes = TRUE, log= T, prior.count = 1)
# # write.csv(normalized.output, file = "03_normalized_data/normalized_output_matrix_log2.csv")
# # 
# # # output object
# # save.image(file = "02_input_data/sfon_wgcna_01_output.RData") # save out existing data 


#### 5. Visualize data ####
#plot using sample IDs
plotMDS(x = my.counts, cex= 0.8) # note that this is supposed to be run on whatever you wrote calcNormFact() to

plotMDS(x = my.counts, cex = 0.8
        , labels = 
#          round(
            interp$date.extract[match(row.names(my.counts$samples), interp$file.name)]
#            )
          )


# #plot using sex
# plotMDS(x = my.counts, cex= 0.8
#         , labels = interp$sex[match(my.counts$samples$files, interp$file.name)])
# #plot using maturity and sex
# plotMDS(x = my.counts, cex= 0.8
#         , labels = paste(
#           interp$sex[match(my.counts$samples$files, interp$file.name)]
#             , interp$poids.sachet.foie[match(my.counts$samples$files, interp$file.name)]
#           , sep = ""))
# 
# # note, this is how matching works:
# interp$sex[match(my.counts$samples$files, interp$file.name)] # matches order 
# interp$sex #see not the same


#### Differential Expression ####
fit <- glmFit(y = my.counts, design = designMat)
lrt <- glmLRT(fit)
edgeR_result <- topTags(lrt)

#?decideTests
deGenes <- decideTestsDGE(lrt, p=0.001)
deGenes <- rownames(lrt)[as.logical(deGenes)]
plotSmear(lrt, de.tags = deGenes)
abline(h=c(-1,1), col =2)
abline(v=cpm.filt, col =2)

deGenes

#### Export Results ####
output <- topTags(lrt, n=nrow(lrt), sort.by="logFC")[[1]]
write.csv(x = output, file = "output.csv")

# Different way:
# group <- colnames(interp$file.name)
# et <- exactTest(my.counts,  = sampleType)


# FINAL FIXES
# need to incorporate the seasonality factor into the DE model