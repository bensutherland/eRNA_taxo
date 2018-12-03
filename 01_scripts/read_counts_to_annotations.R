# Connect read counts to the annotations and plot
# Input: RData from normalization step (post-filter, CPM) and MEGAN (taxonomy ID per amplicon)

# rm(list=ls())

# Choose dataset
datatype <- "eRNA_taxo"

#### 0. Setup (no changes reqd) ####
# Install Packages
#install.packages("RColorBrewer")
library("RColorBrewer")

# Set working directory depending on the dataset
working.dir <- paste("~/Documents/02_eRNA_taxo/", datatype, sep = "")
setwd(working.dir)

## Create a filenames list that contains file names for each dataset (1st obitab output; 2nd MEGAN output)
filenames.list <- list()

# # Root taxonomy level
# filenames.list[["eRNA_taxo"]] <- setNames(object = c(
#                   "08_gx_levels/normalized_output_linear.csv" # counts file
#                 , "16_libs_contig_expr_only_taxonomy_hits-ex.txt") # annot file (TO BE MOVED)
#                                        , nm = c("count", "annot"))

# Family-level
# filenames.list[["eRNA_taxo"]] <- setNames(object = c(
#                 "08_gx_levels/normalized_output_linear.csv" # counts file
#                 , "16_libs_contig_expr_only_taxonomy_hits-ex_Family.txt") # annot file (TO BE MOVED)
#                                       , nm = c("count", "annot"))

# Family-level (and all)
filenames.list[["eRNA_taxo"]] <- setNames(object = c(
  "08_gx_levels/normalized_output_linear.csv" # counts file
  , "16_libs_contig_expr_only_taxonomy_hits-ex_Family_and_all.txt") # annot file (TO BE MOVED)
  , nm = c("count", "annot"))



filenames.list[["eRNA_taxo"]][2]

#### 1.0 Import input data and merge #####
paste("You are analyzing ", datatype, sep = "")
# counts <- read.delim2(file = paste(filenames.list[[datatype]][1], sep = ""), sep = ",") 
load("08_gx_levels/normalized.RData")
str(normalized.output.linear) # this is the counts data, we need to change it though
normalized.output.linear.df <- as.data.frame(normalized.output.linear) # make it a data frame so that we can incorporate the contig ID as a character column
str(normalized.output.linear.df)
normalized.output.linear.df$id <- rownames(normalized.output.linear.df)
str(normalized.output.linear.df)
counts <- normalized.output.linear.df
# now it is ready to continue

# annot <- read.delim2(paste(filenames.list[[datatype]][2], sep = ""), header = F
#                      , col.names = c("id","taxon"))

annot <- read.delim2(filenames.list[["eRNA_taxo"]][2], header = F, col.names = c("id", "taxon")
                     , colClasses = c("character", "character"))
str(annot)


#  ### FRUSTRATING IF NEEDED ###
# # Save the
# rownames(counts) <- counts$X
# counts <- counts[, -1]
# 
# head(counts)
# counts.df <- as.data.frame(counts)
# id <- rownames(counts.df)
# library(dplyr)
# counts <- mutate_all(counts.df, function(x) as.numeric(as.character(x)))
# id <- rownames(counts)
# str(id)
# 
# test <- cbind(counts)
# 
# str(test)
#### END FRUSTRATING IF NEEDED ###

str(counts)
str(annot)
dim(counts)
dim(annot)
head(counts)
head(annot)
names(counts) 
names(annot)

### not needed ? ##
# # This assumes the first column of the counts file is the contig name
# colnames(counts)[1] <- "id"
### end not needed ? ##

# Sort data
counts.sorted <- counts[order(counts$id),]
annot.sorted <- annot[order(annot$id),]

# Merge
data <- merge(x = annot.sorted, y = counts.sorted, by = "id")
names(data)
head(data)
str(data)

# total number of reads dealing with
# Note the difference between normalized data size and non-normalized data size
sum(colSums(data[,3:ncol(data)]))
# per sample
colSums(data[,3:ncol(data)]) # confirm, because this seems greater than the lib sizes earlier..

sum(my.counts$counts[,1]) # see, this is the number of actual reads for first sample, which matches the 'lib.size' within the DGElist. So the normalization must increase this value
sum(my.counts$samples$lib.size) # this is all of the lib sizes


data.df <- data


#### REMOVE UNKNOWNS #####
# Set species to remove (e.g. humans)
remove.from.all <- c("Not assigned", "No hits")
species.remove <- list()
species.remove[["eRNA_taxo"]] <- c(remove.from.all)
species.remove <- species.remove[[datatype]] # Use datatype for removing species
species.remove

# Remove species from data.df
dim(data.df)
data.df <- data.df[ ! data.df$taxon %in% species.remove, ]
dim(data.df) # see how the number of taxa is reduced
### for eRNA_taxo, a ton of data is unassigned or unknown
# goes from 32866 rows (contigs) to 2490!
str(data.df)
sum(data.df[,3:ncol(data.df)]) # still contains 1,891,674 CPM though..

#### 2. Get proportional and count data by taxon per site ####
# Set nulls
sample.tot <- NULL ; sample <- NULL; result <- NULL; result.prop <- NULL ; count.list <- NULL
prop.list <- list(); agg.counts.list <- list()

# Loop to get count by species by site (count.list) and proportion of species by site (prop.list)
for(col in 3:ncol(data.df)) {
  
  # name of the sample this iteration
  sample <- names(data.df[col]) 
  # total number reads in this sample
  sample.tot <- sum(data.df[,col])
  
  # Add 0.0001 to avoid 0 value if sample.tot is 0 (e.g. for controls)
  if(sample.tot==0){
    sample.tot <- 0.00001
  }
  
  # Per sample, aggregate counts by taxon
  result <- aggregate(x = data.df[,col], by = list(data.df$taxon), FUN = sum, na.rm = T)
  
  # Make result proportional by dividing by the amount for that species (2nd column in 'result') 
  result.prop <- (result[,2]/sample.tot)*100 
  
  # Save count values into a list
  count.list[[sample]] <- setNames(result[,2], result[,1]) # pull the value into count.list with names as the first column of result
  
  # Save proportions into a list
  prop.list[[sample]] <- setNames(result.prop, result[,1])
}
  
str(prop.list)
str(count.list)


#### 3. Pull out relevant information from proportional and count data #####
# Proportion data
prop.df <- NULL 

for(i in 1:length(prop.list)){ 
  prop.df <- cbind(prop.df, prop.list[[i]])}

head(prop.df)

# Count data
counts.df <- NULL

for(i in 1:length(count.list)){ 
  counts.df <- cbind(counts.df, count.list[[i]])}

head(counts.df)

##### ADDING SAMPLE NAMES #####
# would be a good place to rename certain columns if necessary
# Incorporate location names
# instead of the following ### NOT WORKING ###
# site.names <- sample.locations[1:length(prop.df[1,])] # Assumes is in same order for the two major types individually (SOG or C3)
site.names <- colnames(data.df[3:ncol(data.df)])
colnames(prop.df) <- site.names # name prop.df
colnames(counts.df) <- site.names # name counts.df
head(prop.df)
head(counts.df)

# Set filenames for saving out
count.output.csv.filename <- paste("09_results/", datatype, "_count_by_taxa.csv", sep = "")
prop.output.csv.filename <- paste("09_results/", datatype, "_prop_by_taxa.csv", sep = "")
# write.csv(x = counts.df, file = count.output.csv.filename)
# write.csv(x = prop.df, file = prop.output.csv.filename)


# Find total numbers of reads mapping (###CPM HERE##)
colnames(counts.df)
sample.reads <- colSums(x = counts.df[, c(1:ncol(counts.df))])


# Filter counts table to remove very low values
min.count <- 10
counts.filtered.df <- counts.df

head(counts.filtered.df)
counts.filtered.df[which(counts.filtered.df < min.count)]

counts.filtered.df[which(counts.filtered.df < min.count)] <- 0
head(counts.filtered.df)

counts.filtered.filename <- paste("09_results/", datatype, "_count_by_taxa_filt_at_", min.count, ".csv", sep = "")
# write.csv(x = counts.filtered.df, file = counts.filtered.filename)

##### 4.0 Prepare plotting (colors) ####
# Prepare palette
#display.brewer.all() # see color options
cols <- brewer.pal(n = 9, name = "Set1")
cols2 <- brewer.pal(n = 8, name = "Set2")
cols3 <- brewer.pal(n = 10, name = "Set3")
cols4 <- brewer.pal(n = 8, name = "Pastel2")
cols5 <- brewer.pal(n = 9, name = "Pastel1")
cols6 <- brewer.pal(n = 11, name = "BrBG")
cols7 <- brewer.pal(n = 10, name = "Paired")
cols8 <- brewer.pal(n = 11, name = "Spectral")
cols9 <- brewer.pal(n = 9, name = "YlOrRd")
cols10 <- brewer.pal(n = 9, name = "YlGnBu")
cols11 <- brewer.pal(n = 9, name = "YlGn")
cols12 <- brewer.pal(n = 9, name = "RdPu")
cols13 <- brewer.pal(n = 9, name = "Purples")
cols14 <- brewer.pal(n = 9, name = "PuRd")
cols15 <- brewer.pal(n = 9, name = "Greys")
cols16 <- brewer.pal(n = 11, name = "RdGy")

palette <- c(cols,cols2,cols3,cols4,cols5,cols6,cols7,cols8,cols9,cols10,cols11,cols12,cols13,cols14,cols15,cols16)
length(palette)

# Randomly select from palette
set.seed(100)
index <- sample(1:nrow(counts.df))
index

# In case need numerous sets
palette.numerous<- rep(x = palette, times = 4)

# Set up in case the number is greater than a single palette
if(length(index) > length(palette)){
  this.palette <- palette.numerous[index]
} else {
  this.palette <- palette[index]
}


# ##### 4.2 Create Legend ####
# # Prepare legend size 
# legend.cex <- c(1, 1, 1, 1, 1, 1) ; names(legend.cex) <- c("C3_16s","C3_COI", "SOG_16s", "C3_val", "SOG_val", "SOG_COI")

# Create dataframe with the taxon and the color
color.index <- cbind(rownames(prop.df), this.palette)
colnames(color.index) <- c("taxon","color")
head(color.index)
color.index.df <- as.data.frame(color.index)
# note this will not work until you move the color codes to character

# Identify which taxa are high proportion in any one sample
min.proport <- 5

# Set null
high.presence.taxa <- NULL

for(i in 1:nrow(prop.df)){
  high.presence.taxa.add <- if(max(prop.df[i,]) > min.proport) { print(rownames(prop.df)[i])}
  high.presence.taxa <- c(high.presence.taxa, high.presence.taxa.add)
}

high.presence.taxa

# Select the rows of the color index for only those high.presence taxa
legend.info <- color.index.df[color.index.df$taxon %in% high.presence.taxa, ]


#### 5. Plot  ####
# filename <- paste("06_output_figures/", datatype, "_read_count_and_prop_by_loc.pdf", sep = "")
# 
# # if want to work interactively, comment out the following line
# pdf(file = filename, width = 10, height = 8)
par(mfrow=c(3,1), mar= c(2.5,4.5,2,1) + 0.2, mgp = c(3.75, 0.75, 0))

# Plot count data
position.info <- barplot(as.matrix(counts.df)
                         , col = this.palette, las = 2, xaxt = "n"
                         , xlim = c(0, ncol(prop.df)+4)
                         , ylab = "Reads")
# axis(side = 1, at = position.info, labels = sample.locations, las = 3, cex.axis = 0.9)

# Plot proportion data
position.info <- barplot(as.matrix(prop.df), col = this.palette
        , xlim = c(0, ncol(prop.df)+4)
        , las = 1
        , ylab = "Proportion (%)"
        , xaxt = "n")

axis(side = 1, at = position.info, 
     labels = site.names, las = 3
     )

# Add information about read counts per sample
# mtext(x = position.info, text = sample.reads
#       , side=3, at = position.info, cex = 0.7)


# Plot Legend, first make blank second plot
plot(1, type = "n", axes = F, xlab = "", ylab = "")

# fix legend info to character text
legend(x = "center", y = "center", legend = legend.info$taxon
        , fill = as.character(legend.info$color)
       #, cex = legend.cex[datatype]
        , ncol = 4)

#dev.off()
#
# Save out as 10 x 8 in portrait

