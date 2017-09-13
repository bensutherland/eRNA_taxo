## Normalization of count data from eXpress output (in table)

#rm(list=ls())

## Install Packages
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
##cdhit95_bbmap: total retained (TRUE) = 20488 # new 20363

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


#### 4. Visualize data ####
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

# Clean up mds plot
sub(pattern = "_.*", replacement = "", x= row.names(my.counts$samples))

plotMDS(x = my.counts, cex = 1.2
        , labels = 
          #          round(
          paste(
            sub(pattern = "_.*", replacement = "", x=
                  interp$file.name[match(row.names(my.counts$samples), interp$file.name)])
          ,  interp$date[match(row.names(my.counts$samples), interp$file.name)]  
          , sep = "_"
          )
        
          
        #            )
)





#### 5. Differential Expression ####
fit <- glmFit(y = my.counts, design = designMat)

lrt <- glmLRT(glmfit = fit, coef = 3)

# Find DE genes for each contrast
lrt.coef1 <- glmLRT(fit, coef = 1) # intercept (all genes)
lrt.coef2 <- glmLRT(fit, coef = 2) # pCO2
lrt.coef3 <- glmLRT(fit, coef = 3) # season
lrt.coef4 <- glmLRT(fit, coef = 4) # pCO2 x season (effect of pCO2 depends on the season)

# Set lrt of choice
LOC <- lrt.coef3
LOC <- lrt.coef2
dim(topTags(LOC, p.value=0.05, n = 15000)) # how many genes DE w/ adj. p-val < 0.05

# # Obtain result
de.result <- topTags(lrt.coef2, p.value=0.05, n = 15000)
dim(de.result)[[1]]
de.result.df <- de.result[[1]]
str(de.result.df)
contig <- rownames(de.result.df)
de.result.out <- cbind(contig, de.result.df)
head(de.result.out)

### 6. Incorporate Annotations ####
# expr.annot <- read.table(file = "../../go_enrichment_annotated_2017-09-11/contig_and_accession.txt", sep = "\t", header = T)
# head(expr.annot)

# Import the sequence annotation separator go_enrichment
expr.annot.all <- read_tsv(file = "../../go_enrichment_annotated_2017-09-11/sequence_annotation.csv", na = "NA")
dim(expr.annot.all)
# convert to df for ease of use
expr.annot.all.df <- as.data.frame(expr.annot.all)
dim(expr.annot.all.df)
colnames(expr.annot.all.df)

# Annotate the full expressed genes
# provides the order that the contigs are in in the dgelist
dgelist.order <- dimnames(my.counts.test[[1]])[1]
dgelist.order <- as.data.frame(dgelist.order[[1]])
str(dgelist.order)

dim(expr.annot.all.df)
colnames(expr.annot.all.df)

head(expr.annot.all.df)[,1]
head(dgelist.order)

str(expr.annot.all.df[,1])
str(dgelist.order)

# Make order a character to match the dataframe it will be merged with
dgelist.order[,1] <- sapply(dgelist.order[,1], as.character)
str(dgelist.order)
colnames(dgelist.order) <- "Contig"

# Confirm data structure type
str(dgelist.order)
str(expr.annot.all.df$Name)

# Order dataframe by target
ordered.annot.df <- expr.annot.all.df[match(dgelist.order$Contig, expr.annot.all.df$Name),]

head(ordered.annot.df$Name)
head(dgelist.order)

dim(ordered.annot.df)
names(ordered.annot.df)


# Export only the expressed genes annotation
write.table(x = ordered.annot.df[,1:2], file = "background_ids.txt", quote = F, sep = "\t", row.names = F, col.names = F)

## Export the background list
# str(dimnames(lrt))
# 
# expr.contigs <- dimnames(lrt)[[1]]
# write.table(x = expr.contigs, file = "all_ids.txt" , quote = F, sep = "\t"
#             , row.names = F, col.names = FALSE)


annotated.de <- merge(de.result.out, expr.annot.all.df, by.x = "contig", by.y = "Name" )
dim(de.result.out)
dim(annotated.de)
str(annotated.de)

# write annotated de
write.table(x = annotated.de, file = "annotated_de.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = annotated.de, file = "annotated_de_season.txt", quote = F, sep = "\t", row.names = F, col.names = T)

#### 7. Export Results ####
write.table(x = x, file = "differential_genes.txt", quote = F, sep = "\t"
            , col.names = TRUE, row.names = F)
write.table(x = contig, file = "significant_ids.txt", quote = F, sep = "\t", row.names = F, col.names = F)

# Export only expressed genes:
expr.contigs <- dimnames(lrt)[[1]]
write.table(x = expr.contigs, file = "all_ids.txt" , quote = F, sep = "\t"
          , row.names = F, col.names = FALSE)


# Temp (working)
# to plot a single gene:
boxplot(my.counts[[1]]["contig-60_2892_length_1823_read_count_49",] ~ interp$Range)

# to view contrasts
lrt.coef2
