# Assess the assembly results and plot
# Must have no inputlibrary size for the oddball sample (the one that was assembled individually then collapsed in)

setwd("~/Documents/02_eRNA_oyster")

assembly.records <- read.csv("eRNA_metatxome_assembly_records_2017-10-11.csv")

names(assembly.records)

par(mfrow=c(2,1), mar= c(4,3,0.5,1) + 0.2, mgp = c(2,0.75,0))

#### Plot Assembly Stats ####
ymax <- max(assembly.records$num.ctg.out/1000000, assembly.records$tot.leng.ctg.out/1000000000, na.rm=T)
xmax <- max(assembly.records$tot.Gb.in, na.rm = T)

odd.lib.x <- 47

plot(x = assembly.records$tot.Gb.in
     , xlim = c(0, xmax + xmax * 0.1)
     , ylim = c(0, ymax + ymax * 0.1)
     , type = "n", las = 1
     , xlab = "Input data size (Gb)"
     , ylab = "")

# plot million output contigs
million.ctgs <- assembly.records$num.ctg.out/1000000
points(x = assembly.records$tot.Gb.in, y = million.ctgs , col = "red", pch = 15)
lines(x = assembly.records$tot.Gb.in, y = million.ctgs, col = "red", lty=3)

# include eight lib joined
points(x = odd.lib.x, y = million.ctgs[is.na(assembly.records$tot.reads.in)], col = "red", pch = 16)

# plot Gb total length
tot.len.Gb <- assembly.records$tot.leng.ctg.out/1000000000
points(x = assembly.records$tot.Gb.in, y =  tot.len.Gb, col = "black", pch = 15)
lines(x = assembly.records$tot.Gb.in, y = tot.len.Gb, col = "black", lty=3)

#include eight lib joined
points(x = odd.lib.x, y = tot.len.Gb[is.na(assembly.records$tot.reads.in)], col = "black", pch = 16)

legend(x = "bottomright", legend = c("# ctg (M)", "tot len (Gb)"), fill = c("red","black"))

# Add number libraries
text.level <- 8.4
num.libs <- c(1,2,4,8,12,16)
text(x = assembly.records$tot.Gb.in[0:length(num.libs)], y = text.level, labels = num.libs)
text(x = 55, y = text.level, labels = "# Libs")
text(x = odd.lib.x, y = text.level, labels = 8)

abline(v = 44, lty = 2)

#### Plot percentages ####
names(assembly.records)


# Plot
ymax <- 100
xmax <- max(assembly.records$tot.Gb.in, na.rm = T)

plot(x = assembly.records$tot.Gb.in
     , xlim = c(0, xmax + xmax * 0.1)
     , ylim = c(0, ymax + ymax * 0.1)
     , type = "n", las = 1
     , xlab = "Input data size (Gb)"
     , ylab = "Percent (%)")

points(x = assembly.records$tot.Gb.in, y = assembly.records$s15.align.conc.0, col = "darkgreen", pch = 15)
lines(x = assembly.records$tot.Gb.in, y = assembly.records$s15.align.conc.0, col = "darkgreen", lty=3)
points(x = assembly.records$tot.Gb.in, y = assembly.records$s15.align.conc.1, col = "blue", pch = 15)
lines(x = assembly.records$tot.Gb.in, y = assembly.records$s15.align.conc.1, col = "blue", lty=3)
points(x = assembly.records$tot.Gb.in, y = assembly.records$s15.align.conc.morethan1, col = "darkgrey", pch = 15)
lines(x = assembly.records$tot.Gb.in, y = assembly.records$s15.align.conc.morethan1, col = "darkgrey", lty=3)

# include eight lib joined
points(x = odd.lib.x, y = assembly.records$s15.align.conc.0[is.na(assembly.records$tot.reads.in)], col = "darkgreen", pch = 16)
points(x = odd.lib.x, y = assembly.records$s15.align.conc.1[is.na(assembly.records$tot.reads.in)], col = "blue", pch = 16)
points(x = odd.lib.x, y = assembly.records$s15.align.conc.morethan1[is.na(assembly.records$tot.reads.in)], col = "darkgrey", pch = 16)

abline(v = 44, lty = 2)

legend(x = "bottomright", legend = c("Not align", "Align 1", "Align > 1"), fill = c("darkgreen","blue", "darkgrey"))


# save out as 8 x 5.5 in portrait pdf

