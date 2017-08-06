### GRAPHICALLY REPRESENT SILHOUETTE SCORES, 
# one input parameter - prefix of silhouette output files
library(gridExtra)
library(RColorBrewer)

arg <-commandArgs()
ind <- read.table(paste(arg[6],"_ind.sil",sep=""), header = F)
pop <- read.table(paste(arg[6],"_pop.sil",sep=""), header = F)
colnames(pop)=c("pop","mean_SIL","nemesis","nemesis_freq")

palette <-c(brewer.pal(11,"BrBG"), brewer.pal(11, "PRGn"), brewer.pal(11, "RdGy"))


cols <- c()

for (i in 1:length(table(ind$V1)))
{
	cols <- c(cols, rep(palette[i],table(ind$V1)[i]))
}
pdf(paste(arg[6],"_ind.pdf",sep=""))
barplot(ind$V2, col = cols, border = NA)

legend("topright", legend = names(table(ind$V1)), cex = 0.5, col = palette, pch = 15,  box.lwd = 0)
dev.off()


pdf(paste(arg[6],"_pop.pdf",sep=""))
barplot(pop$mean_SIL, col = palette, border = NA)

legend("topright", legend = pop$pop, cex = 0.5, col = palette, pch = 15,  box.lwd = 0)
plot(NA, ylab = "", xlab = "", ylim = c(0,1), xlim = c(0,1),xaxt = 'n', yaxt = 'n', frame = F)
grid.table(pop)
dev.off()
