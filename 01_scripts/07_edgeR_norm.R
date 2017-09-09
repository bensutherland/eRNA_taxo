## Normalization of count data from eXpress output (in table)

#rm(list=ls())

## Install Packages
#source("http://bioconductor.org/biocLite.R")
#biocLite("edgeR")
require("edgeR")

#biocLite("locfit")
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
dim(my.counts) 
##bbmap_only: total unique tags =    4528162
##cdhit95_bbmap: total unique tags = 4081958 

## Backup stuff
# my.counts.bk <- my.counts # just in case
# my.counts <- my.counts.bk # CAUTION go backwarks

# Set up DGEList
rownames(my.counts) <- my.counts[,1]
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
##bbmap_only: total retained (TRUE) =    
##cdhit95_bbmap: total retained (TRUE) = 20488 

# subset DGEList
my.counts <- my.counts[keep, , keep.lib.sizes=FALSE] #keep.lib.sizes = T retains original lib sizes, otherwise recomputes w remaining tags
dim(my.counts)

#### 3. Normalization ####
# Use TMM normalization, as it takes into account highly expressed genes that may take up sequencing rxn and make other genes look down-reg.
my.counts <- calcNormFactors(my.counts, method = c("TMM"))
my.counts$samples

# Plot norm.factors by library size
plot(my.counts$samples$norm.factors ~ my.counts$samples$lib.size)

#### 3. Create design matrix ####
# Use approximated variables in bins
binary.pCO2 <- interp$Range 
binary.season <- interp$season

#missing s13
## binary.pCO2 <- interp$Range[interp$file.name!="S13"] # currently missing one
## binary.season <- interp$season[interp$file.name!="S13"] # currently missing one

# Build a design matrix
#designMat <- model.matrix(~binary.pCO2)
designMat <- model.matrix(~binary.pCO2 * binary.season)
#designMat <- model.matrix(~binary.season)

# Estimate dispersions (measure inter-library variation per tag)
# IT appears that estimateDisp is for simpler models, whereas estimateGLMCommon etc. is for when doing glms
my.counts <- estimateDisp(my.counts, design=designMat) # note that this can use a design matrix when provided 
#"qCML method is only applicable on datasets with a single factor design
#since it fails to take into account the effects from multiple factors in a more complicated
#experiment."


# my.counts <- estimateGLMCommonDisp(my.counts, design=designMat)
# my.counts <- estimateGLMTrendedDisp(my.counts, design=designMat)
# my.counts <- estimateGLMTagwiseDisp(my.counts, design=designMat)

summary(my.counts$prior.df) # est. overall var. across genome for dataset
sqrt(my.counts$common.disp) #coeff of var, for biol. var

plotBCV(my.counts)


#### 5. Visualize data ####
#plot using sample IDs
plotMDS(x = my.counts, cex= 0.8) # note that this is supposed to be run on whatever you wrote calcNormFact() to

# Note that it is the date of the sample that explains PC1
plotMDS(x = my.counts, cex = 0.8
        , labels = 
#          round(
            interp$date[match(row.names(my.counts$samples), interp$file.name)]
#            )
          )

# # note, this is how matching works:
# interp$sex[match(my.counts$samples$files, interp$file.name)] # matches order 
# interp$sex #see not the same

#### 6. Differential Expression ####
fit <- glmFit(y = my.counts, design = designMat)

lrt <- glmLRT(glmfit = fit, coef = 3)

# Find DE genes for each contrast
lrt.coef1 <- glmLRT(fit, coef = 1) # intercept (all genes)
lrt.coef2 <- glmLRT(fit, coef = 2) # pCO2
lrt.coef3 <- glmLRT(fit, coef = 3) # season
lrt.coef4 <- glmLRT(fit, coef = 4) # pCO2 x season (effect of pCO2 depends on the season)

# Set lrt of choice
LOC <- lrt.coef4
dim(topTags(LOC, p.value=0.05, n = 15000)) # how many genes DE w/ adj. p-val < 0.05

# # Obtain result
# edgeR_result <- topTags(lrt.coef2, p.value=1, n = 15000)
# dim(edgeR_result)

# ## Alternate method, usign decideTests
# #?decideTests
# deGenes <- decideTestsDGE(lrt, p=0.001)
# deGenes <- rownames(lrt)[as.logical(deGenes)]
# plotSmear(lrt, de.tags = deGenes)
# abline(h=c(-1,1), col =2)
# abline(v=cpm.filt, col =2)
# 
# deGenes

#### Export Results ####
output <- topTags(lrt, n=nrow(lrt), sort.by="logFC")[[1]]
write.csv(x = output, file = "output.csv")

# Different way:
# group <- colnames(interp$file.name)
# et <- exactTest(my.counts,  = sampleType)


# FINAL FIXES
# need to incorporate the seasonality factor into the DE model

# Export only expressed genes:
expr.contigs <- dimnames(lrt)[[1]]
write.table(x = expr.contigs, file = "expr_contigs.csv" , quote = F, sep = ","
          , row.names = F, col.names = FALSE)

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
