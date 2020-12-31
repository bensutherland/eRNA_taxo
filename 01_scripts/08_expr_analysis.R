# Differential Expression Analysis and Data Plotting
# 2020-12-31
# B. Sutherland

# rm(list=ls())

# Set working directory
setwd("/hdd/20_other_research/metatranscriptomics/eRNA_taxo")

# Load data from normalization script (i.e., 01_scripts/07_edgeR_norm.R)
load(file = "08_gx_levels/normalized.RData")

### 1. Bring in annotation and make sure matches gene expr object ####
# note: will result in the expressed genes having annnotation associated

# Import sequence annotation
expr.annot.all.df <- read.delim2(file = "06_metatranscriptome/sequence_annotation.txt", header = T, sep = "\t")
expr.annot.all.df <- as.data.frame(expr.annot.all.df, stringsAsFactors = FALSE)
colnames(expr.annot.all.df)
expr.annot.all.df$Name[1:6]
dim(expr.annot.all.df)
length(unique(expr.annot.all.df$Name))

# In sequence annotation, retain contig names only up to the first space to match with gene expr obj.
expr.annot.all.df$Name <- gsub(pattern = "\ .*", replacement = "", x = expr.annot.all.df$Name)
length(unique(expr.annot.all.df$Name)) # should have remained the same length as prior to removing everything after the first space

# Identify contig order in the gene expression obj (i.e., DGElist), save as a df
dgelist.order <- dimnames(my.counts[[1]])[1]
dgelist.order <- as.data.frame(x = dgelist.order[[1]], stringsAsFactors = FALSE)
colnames(dgelist.order) <- "Contig"
str(dgelist.order) # confirm is character
head(dgelist.order)

# Order dataframe by target
ordered.annot.df <- expr.annot.all.df[match(dgelist.order$Contig, expr.annot.all.df$Name), ]
head(ordered.annot.df$Name) #to confirm
head(dgelist.order)

# View annotation df
dim(ordered.annot.df)
names(ordered.annot.df)
str(ordered.annot.df)

# convert to character
ordered.annot.df <- data.frame(lapply(ordered.annot.df, as.character), stringsAsFactors =F)
str(ordered.annot.df)
head(ordered.annot.df)

# Export the annotation
write.table(x = ordered.annot.df, file = "09_results/background_ids.txt", quote = F, sep = "\t", row.names = F, col.names = T)


#### 2. Visualize PCA/MDS ####
# Explore data and colour by selected variable
colnames(interp) # variable options
variable <- "temp" # choose variable
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
pdf(file = "09_results/PCA_PCO2_date_season.pdf", width = 12.5, height = 8.5)
trait <- "season"
MDS.plot.res <- plotMDS(x = my.counts, cex = 0.8, gene.selection = "common"
                        , col = as.numeric(interp[match(row.names(my.counts$samples), interp$file.name), trait])
                        , labels = 
                          paste(sep = "_" 
                                # sample name, remove everything after first underscore
                                , sub(pattern = "_.*", replacement = "", x=interp$file.name[match(row.names(my.counts$samples), interp$file.name)])
                                
                                # pCO2 level
                                , round(interp[match(row.names(my.counts$samples), interp$file.name), "pCO2"], digits = 1)
                                
                                # date
                                , interp[match(row.names(my.counts$samples), interp$file.name), "date"]
                          )
)

legend(x="bottomright", legend=c("early summer", "late summer/fall"), fill = c("black","red"), cex = 0.8)
#text(x = 1.5, y = 2, labels = "Sample_pCO2", cex = 1.2)
# save out as 6.5 x 6.5
dev.off()


# ### Alternative mapping method with symbols
# plotMDS(x = my.counts, cex = 0.8
#         , pch = as.numeric(interp[match(row.names(my.counts$samples), interp$file.name), trait])
#         )

# ### Alternative retain coordinates to use later
# plot.data <- plotMDS(x = my.counts)
# then use $x and $y to place text


#### 3. Differential expression analysis ####
# Use approximated variables in bins
binary.pCO2 <- interp$Range
binary.season <- interp$season

# Build a design matrix
designMat <- model.matrix(~binary.pCO2 * binary.season)

# Estimate dispersions (measures inter-library variation per tag)
# estimateDisp for simpler models; estimateGLMCommon (etc) is for when doing glms
my.counts <- estimateGLMCommonDisp(my.counts, design=designMat)
my.counts <- estimateGLMTrendedDisp(my.counts, design=designMat)
my.counts <- estimateGLMTagwiseDisp(my.counts, design=designMat)

summary(my.counts$prior.df) # est. overall var. across genome for dataset
sqrt(my.counts$common.disp) # coeff of var, for biol. var
plotBCV(my.counts)

# create glm fit
fit <- glmFit(y = my.counts, design = designMat)
#lrt <- glmLRT(glmfit = fit, coef = 3)

# Here is the design matrix to identify which intercept you are interested in
head(designMat)

# Find DE genes for each contrast
lrt.coef1 <- glmLRT(fit, coef = 1) # intercept (all genes)
lrt.coef2 <- glmLRT(fit, coef = 2) # pCO2
lrt.coef3 <- glmLRT(fit, coef = 3) # season
lrt.coef4 <- glmLRT(fit, coef = 4) # pCO2 x season (effect of pCO2 depends on the season)

# Obtain results
de.result.pCO2 <- topTags(lrt.coef2, p.value=0.05, adjust.method = "BH", n = 15000) #pCO2
de.result.season <- topTags(lrt.coef3, p.value=0.05, adjust.method = "BH", n = 15000) #season
de.result.ifx <- topTags(lrt.coef4, p.value=0.05, adjust.method = "BH", n = 15000) #season

nrow(de.result.pCO2) # 720
nrow(de.result.season) # 2765
nrow(de.result.ifx) # 0

head(de.result.pCO2)
head(de.result.season)

# Put results into a list of lists
de_results.list <- list()
de_results.list[["pCO2"]] <- de.result.pCO2
de_results.list[["season"]] <- de.result.season
de_results.list[["ifx"]] <- de.result.ifx

str(de_results.list)

# Format the items in the list
result.oi <- NULL; result.oi.df <- NULL; all_results.df <- NULL

for(d in 1:length(de_results.list)){
  
  # Ignore any contrasts without significant results
  if(length(de_results.list[[d]])!=0){
  
    result.oi <- de_results.list[[d]]
    
    # Pull the df out of the list slot
    result.oi.df <- result.oi[["table"]]
    
    # Bring vector of contigs into the table
    result.oi.df$contig <- rownames(result.oi.df)
    
    # Bring contrast type into the df
    result.oi.df$contrast <- rep(x = result.oi$comparison, times = nrow(result.oi.df))
    
    # Reorder
    result.oi.df <- result.oi.df[,c("contig", "contrast", "logFC", "logCPM", "LR", "PValue", "FDR")]
    
    all_results.df <- rbind(all_results.df, result.oi.df)
    
    # Reporting
    print(paste0("The contrast *", names(de_results.list)[d], "* has been added to the output. "))
  
  }else{
    
    print(paste0("The contrast *", names(de_results.list)[d], "* has no results and will be discarded. "))
    
  }
  
}

str(all_results.df) # contig should be character

# Merge results and annotation
annotated.de <- merge(x = all_results.df, y = ordered.annot.df, by.x = "contig", by.y = "Name", all.x = TRUE)
dim(all_results.df)
dim(annotated.de)
str(annotated.de)

# Sort by contrast type
annotated.de <- annotated.de[order(annotated.de$contrast, annotated.de$logFC), ]
head(annotated.de)


#### 4. Export Results ####
# Save out the differential gene list
write.table(x = annotated.de, file = "09_results/diff_genes_annotated.txt", quote = F, sep = "\t", row.names = F, col.names = T)



###### EXTRA MATERIALS ######

#### 5. Plotting GOIs #####
GOI1 <- "contig-60_1777509" # Terminase, large subunit
GOI2 <- "contig-60_1307617" # Tail sheath protein (Vibrio phage)
GOI3 <- "contig-60_557998" # DNA-directed RNA polymerase
GOI4 <- "contig-60_495817" # "Sliding-clamp-loader gp44 subunit"

# # Single Gene
# GOI1 <- GOI2
# boxplot(my.counts[[1]][GOI,] ~ interp$season, main = GOI)

gene.names <- c("Terminase, large subunit"
                , "Tail sheath protein"
                , "DNA-directed RNA polymerase"
                , "Sliding-clamp-loader gp44 subunit")

# Plot 2 x 2 individual genes
par(mfrow=c(2,2), mar= c(2,3,2,1) + 0.2, mgp = c(2,0.75,0))


# These should not be extracted from my.counts, but rather from the normalized cpm data for plotting
# extracted from my.counts (incorrect)
# GOIs <- c(GOI1, GOI2, GOI3, GOI4)
# for(goi in 1:length(GOIs)){
#   boxplot(my.counts[[1]][GOIs[goi],] ~ interp$season
#           , main = gene.names[goi]
#           , las = 1
#           , ylab = "read counts")
# }
# save out as 11 x 7 in portrait

# my.counts[[1]][1:5,1:5]
# colSums(my.counts[[1]])

# extracted from cpm-normalized
GOIs <- c(GOI1, GOI2, GOI3, GOI4)
for(goi in 1:length(GOIs)){
  boxplot(normalized.output.linear[GOIs[goi],] ~ interp$season
          , main = gene.names[goi]
          , las = 1
          , ylab = "read counts (CPM)")
}
