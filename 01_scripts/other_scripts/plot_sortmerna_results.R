# Use this to plot rRNA vs mRNA reads
setwd("~/Documents/02_eRNA_oyster")
#options(scipen=999)
options(scipen=0)

read.data <- read.table(file = "overview_sortmerna_results_2017-07-13.csv", sep = ",", header = T)
head(read.data)

non.rRNA.reads <- read.data$Fail_E.val
sample.ids <- read.data$SampleID
total.reads <- read.data$Reads

data <- as.data.frame(cbind(
  #sample.ids, 
  non.rRNA.reads, total.reads ))
row.names(data) <- sample.ids

# Plot
barplot(data$total.reads/1000000, las = 1
        , ylab = "Million Reads", col = "lightgrey"
        , names = sample.ids, cex.names = 0.8
        )
barplot(data$non.rRNA.reads/1000000, add = T, col = "darkgrey"
        , yaxt = "n")
legend(x = "topleft", legend = c("all reads", "non-rRNA reads")
       , fill = c("lightgrey","darkgrey"))
