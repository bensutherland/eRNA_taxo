setwd("~/Documents/02_eRNA_taxo/eRNA_taxo")

alignment.rate <- read.table("07_mapped/overall_alignment_rate.txt")

summary(alignment.rate)

align.concord <- read.table("07_mapped/align_concord_exact_once.txt")

summary(align.concord)

align.concord.more.than.one <- read.table("07_mapped/align_concord_more_than_once.txt")

summary(align.concord.more.than.one)
