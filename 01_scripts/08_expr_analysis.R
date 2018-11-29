# Differential Expression Analysis and Data Plotting

# rm(list=ls())

# Load data from normalization script
load(file = "08_gx_levels/normalized.RData")

### 1. Incorporate annotation ####
# Import sequence annotation
expr.annot.all <- read.delim2(file = "06_metatranscriptome/sequence_annotation.txt", header = T, sep = "\t")
head(expr.annot.all)
head(expr.annot.all$Name) # need to shorten to the first space
expr.annot.all$Name <- gsub(pattern = "\ .*", replacement = "", x = expr.annot.all$Name) # keep only the first part of the name to match
dim(expr.annot.all)
expr.annot.all.df <- as.data.frame(expr.annot.all) # convert to df for ease of use
dim(expr.annot.all.df)
colnames(expr.annot.all.df)
head(expr.annot.all.df)


# Combine the DGElist with the annotation
#  NOT DONE YET 

##### is this section needed?
# provides the order that the contigs are in in the dgelist
dgelist.order <- dimnames(my.counts[[1]])[1]
dgelist.order <- as.data.frame(dgelist.order[[1]])
str(dgelist.order)

dgelist.order[,1] <- sapply(dgelist.order[,1], as.character) # Order into character for df matching
str(dgelist.order)
colnames(dgelist.order) <- "Contig"

str(dgelist.order) # Confirm data structure type
str(expr.annot.all.df$Name)

# Order dataframe by target
ordered.annot.df <- expr.annot.all.df[match(dgelist.order$Contig, expr.annot.all.df$Name),]
head(ordered.annot.df$Name) #to confirm
head(dgelist.order)

dim(ordered.annot.df)
names(ordered.annot.df)
str(ordered.annot.df)

# convert to character
ordered.annot.df <- data.frame(lapply(ordered.annot.df, as.character), stringsAsFactors =F)
str(ordered.annot.df)

# result: ordered.annot.df having annotation of all expressed genes

# Aside, for go_enrichment method of GO enrich
## Export the background list
# str(dimnames(lrt))
#
# expr.contigs <- dimnames(lrt)[[1]]
# write.table(x = expr.contigs, file = "all_ids.txt" , quote = F, sep = "\t"
#             , row.names = F, col.names = FALSE)
##### is this section needed? (end)


#### 2. Visualize PCA/MDS ####
# Exploration plot (with sample ID short and selected variable)
colnames(interp)
variable <- "Calcite.om"
# Note: if variable is numeric, you may need to uncomment out the digits argument below

plotMDS(x = my.counts, cex = 0.8
        , labels = paste(sep = "_",
                         sub(pattern = "_.*", replacement = "", 
                             x=interp$file.name[match(row.names(my.counts$samples), interp$file.name)])

                         ,  
#                         round(
                         interp[match(row.names(my.counts$samples), interp$file.name), variable]  
#                         , digits = 2)
        )
)


# Plot for report (sampleID_pCO2 colored by season)
# Note: using gene.selection = common instead of pairwise means that the same genes are used to distinguish all samples
pdf(file = "09_results/PCA_plot.pdf", width = 12.5, height = 8.5)
MDS.plot.res <- plotMDS(x = my.counts, cex = 0.8, gene.selection = "common"
                        , col = as.numeric(interp[match(row.names(my.counts$samples), interp$file.name), trait])
                        , labels = 
                          paste(sep = "_" 
                                , sub(pattern = "_.*", replacement = ""
                                      , x=interp$file.name[match(row.names(my.counts$samples), interp$file.name)])
                                , round(interp[match(row.names(my.counts$samples), interp$file.name), "pCO2"], digits = 1)
                                , interp[match(row.names(my.counts$samples), interp$file.name), "date"]
                          )
)

legend(x="bottomright", legend=c("early summer", "late summer/fall"), fill = c("black","red"), cex = 0.8)
#text(x = 1.5, y = 2, labels = "Sample_pCO2", cex = 1.2)
# save out as 6.5 x 6.5
dev.off()

### Alternate options:
## e.g. use symbols
plotMDS(x = my.counts, cex = 0.8
        , pch = as.numeric(interp[match(row.names(my.counts$samples), interp$file.name), trait])
        )

## e.g. save coordinates so you can add text later
# plot.data <- plotMDS(x = my.counts)
# then use $x and $y to place text


##### Alternate PCA approach ####
y <- cpm(my.counts, log = T, prior.count = 3) # returns log2 vals
pca <- prcomp( t( y ), scale  = T)

# Use package factoextra for visualization
# install.packages("factoextra")
library("factoextra")
fviz_eig(pca) # view the Percent variance explained (first three matter)
fviz_pca_ind(pca, repel = T, axes = c(1,2))
fviz_pca_ind(pca, repel = T, axes = c(2,3))

# find top loading genes?
res.var <- get_pca_var(pca)
head(res.var$contrib)
boxplot(res.var$contrib[,1]) #for example, boxplot the top loading genes
summary(res.var$contrib[,1])
names(which(res.var$contrib[,1] > 0.01)) # e.g. top loading..


#### 3. Differential expression analysis ####
# Use approximated variables in bins
binary.pCO2 <- interp$Range 
binary.season <- interp$season

# Build a design matrix
designMat <- model.matrix(~binary.pCO2 * binary.season)

# Estimate dispersions (measure inter-library variation per tag)
# to note, then delete #
# estimateDisp for simpler models; estimateGLMCommon (etc) is for when doing glms
my.counts <- estimateGLMCommonDisp(my.counts, design=designMat)
my.counts <- estimateGLMTrendedDisp(my.counts, design=designMat)
my.counts <- estimateGLMTagwiseDisp(my.counts, design=designMat)

summary(my.counts$prior.df) # est. overall var. across genome for dataset
sqrt(my.counts$common.disp) #coeff of var, for biol. var

plotBCV(my.counts)

# create glm fit
fit <- glmFit(y = my.counts, design = designMat)
#lrt <- glmLRT(glmfit = fit, coef = 3)

designMat
# Find DE genes for each contrast
lrt.coef1 <- glmLRT(fit, coef = 1) # intercept (all genes)
lrt.coef2 <- glmLRT(fit, coef = 2) # pCO2
lrt.coef3 <- glmLRT(fit, coef = 3) # season
lrt.coef4 <- glmLRT(fit, coef = 4) # pCO2 x season (effect of pCO2 depends on the season)

# Obtain results
de.result.pCO2 <- topTags(lrt.coef2, p.value=0.05, n = 15000) #pCO2
de.result.season <- topTags(lrt.coef3, p.value=0.05, n = 15000) #season
de.result.ifx <- topTags(lrt.coef4, p.value=0.05, n = 15000) #season

nrow(de.result.pCO2) #747
nrow(de.result.season) #2792
nrow(de.result.ifx) # 3

# Choose result
#de.result <- de.result.pCO2
de.result <- de.result.season
#de.result <- de.result.ifx

# Final formatting
dim(de.result)[[1]]
de.result.df <- de.result[[1]]
str(de.result.df)
contig <- rownames(de.result.df)
de.result.out <- cbind(contig , de.result.df)
de.result.out$contig <- as.character(de.result.out$contig)
head(de.result.out)
str(de.result.out)

# Bring in annotation onto the de.result
annotated.de <- merge(de.result.out, ordered.annot.df, by.x = "contig", by.y = "Name" )
dim(de.result.out)
dim(annotated.de)
str(annotated.de)


#### 4. Export Results ####
# Choose appropriate name and save out the differential gene list
#write.table(x = annotated.de, file = "annotated_de_pCO2.txt", quote = F, sep = "\t", row.names = F, col.names = T)
write.table(x = annotated.de, file = "annotated_de_season.txt", quote = F, sep = "\t", row.names = F, col.names = T)
#write.table(x = annotated.de, file = "annotated_de_ifx.txt", quote = F, sep = "\t", row.names = F, col.names = T)

# For go_enrich version 
# write.table(x = contig, file = "significant_ids.txt", quote = F, sep = "\t", row.names = F, col.names = F)
# # Export only expressed genes:
# expr.contigs <- dimnames(lrt)[[1]]
# write.table(x = expr.contigs, file = "all_ids.txt" , quote = F, sep = "\t"
#           , row.names = F, col.names = FALSE)

# Export all the expressed genes annotation
write.table(x = ordered.annot.df[,1:2], file = "background_ids.txt", quote = F, sep = "\t", row.names = F, col.names = F)


# Aside, to view contrasts
# lrt.coef2



#### HERE BE DRAGONS #####

#### 8. Visualize single genes ####
# to plot a single gene:
boxplot(my.counts[[1]]["contig-60_2892_length_1823_read_count_49",] ~ interp$Range )
boxplot(my.counts[[1]]["contig-60_1051284_length_201_read_count_158",] ~ interp$season)
boxplot(my.counts[[1]]["contig-60_1051284_length_201_read_count_158",] ~ interp$Range)

# example of one DE by both:
boxplot(my.counts[[1]]["contig-60_100194_length_278_read_count_9",] ~ interp$season * interp$Range)

# GOI style:
#GOI <- "contig-60_653007_length_231_read_count_2"
GOI <- "contig-60_10002_length_800_read_count_47"
boxplot(my.counts[[1]][GOI,] ~ interp$season * interp$Range, main = GOI)


# Top candidates:
#GOI1 <- "contig-60_104578_length_274_read_count_12" # DNA-directed DNA polymerase
GOI1 <- "contig-60_11283_length_1039_read_count_71" # Terminase, large subunit #Good #*#
#GOI1 <- "contig-60_1356_length_1536_read_count_48" # Capsid, very low expressed
#GOI1 <- "contig-60_169441_length_335_read_count_12" # Terminase again, good
GOI2 <- "contig-60_950163_length_209_read_count_4" # DNA-directed DNA polymerase #*#, higher expr in high summer
#GOI1 <- "contig-60_65888_length_522_read_count_6" # exonuclease #*#

GOI <- GOI1
boxplot(my.counts[[1]][GOI,] ~ interp$season * interp$Range, main = GOI, las = 1)


# IFX
GOI3 <- "contig-60_15091_length_922_read_count_56" # clamp
#GOI3 <- "contig-60_65888_length_522_read_count_6" # exonuclease probable (same as others)


GOI <- GOI3
boxplot(my.counts[[1]][GOI,] ~ interp$season * interp$Range, main = GOI, las = 1) # v nice.. same #*#

# MFX season
GOI4 <- "contig-60_901732_length_210_read_count_0" # CRISPR-assoc #*#
#GOI4 <- "contig-60_715759_length_223_read_count_4" # capsid
#GOI4 <- "contig-60_118_length_4120_read_count_834" # Tail tubular protein

GOI <- GOI4
boxplot(my.counts[[1]][GOI,] ~ interp$season * interp$Range, main = GOI, las = 1)

gene.names <- c("Terminase, large subunit"
                , "DNA-directed DNA polymerase"
                , "Clamp loader large subunit"
                , "CRISPR-assoc. endonuc. Cas9")

# Plot 2 x 2 individual genes
par(mfrow=c(2,2), mar= c(2,3,2,1) + 0.2, mgp = c(2,0.75,0))

GOIs <- c(GOI1, GOI2, GOI3, GOI4)
for(goi in 1:length(GOIs)){
  boxplot(my.counts[[1]][GOIs[goi],] ~ interp$season * interp$Range
          , main = gene.names[goi]
          , las = 1
          , ylab = "read counts")
}
# save out as 11 x 7 in portrait
