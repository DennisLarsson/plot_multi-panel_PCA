#version 2.11
library(adegenet)
#if you dont have adegenet installed, install using: install.packages("adegenet")

setwd("/home/biogeoanalysis/RAD/spicatumGroup/06populations_50miss_mac2/")
infile<-"spicgrp.stru"
outputfile.name <-"spicgrp_popLevel"
popmap="/home/biogeoanalysis/RAD/spicatumGroup/popmap_spicGroup_sub_Grp"      #give the name of the popmap if it is in the work directory or the full path if it is somewhere else
number.indv=89
number.loci=41912

ld.pos.PC12 = "topright" #where the legend will be placed in the PC1-2 plot: topright, topleft, bottomright, bottomleft
ld.pos.PC13 = "topleft" #where the legend will be placed in the PC1-3 plot: topright, topleft, bottomright, bottomleft

pop <-read.delim(popmap, header = FALSE, as.is=T)
#make sure to edit n.ind to the number of individuals your dataset has and n.loc to the number of loci you have.
object1 <- read.structure(infile, n.ind = number.indv, n.loc = number.loci, onerowperind = FALSE, col.lab = 1, col.pop = 2, row.marknames = 1, 
                          NA.char = "-9", ask = FALSE)

X <- tab(object1, freq = TRUE, NA.method = "mean")

pca1 <- dudi.pca(X, scale = FALSE, scannf = FALSE, nf = 3)

powerVal1 = round((pca1$eig[1]/sum(pca1$eig))*100,digits = 2)
powerVal2 = round((pca1$eig[2]/sum(pca1$eig))*100,digits = 2)
powerVal3 = round((pca1$eig[3]/sum(pca1$eig))*100,digits = 2)
val<-c(powerVal1,powerVal2,powerVal3)


# introduce function which goes through the popmap and indexes each uniq occurence of popname, then add this number in a new list
uniq.pop <- unique(pop[,2])
pop.shape <- data.frame(shape=1:length(uniq.pop),row.names = uniq.pop)
shape=1
i=1
while (i <= length(uniq.pop)) {
  pop.shape[i,1]=shape
  i=i+1
  shape=shape+1
  if (shape > 18) {
    shape=1
  }
}
color.pal <- c("brown", "blue", "red", "green","orange", "purple")
pop.features <- data.frame(color=1:length(pop[,2]),shape=1:length(pop[,2]),row.names = pop[[1]])

for (i in 1:length(pop[,2])) {
  pop.match=match(pop[i,2],uniq.pop)
  pop.features[i,1]=color.pal[pop.match]
  pop.features[i,2]=pop.shape[pop.match,1]
}

plot.PCA <- function(x,y,pos,plot.cex,plot.pt.cex,plot.ncol,prefix) {
  plot(pca1$li[1,x], pca1$li[1,y], col=pop.features[1,1], pch=pop.features[1,2], 
       xlab=paste("PC",x," - ",val[x],"%",sep=""), ylab=paste("PC",y," - ",val[y],"%",sep=""), 
       xlim=c(min(pca1$li[,x]),max(pca1$li[,x])), ylim=c(min(pca1$li[,y]),max(pca1$li[,y])), 
       main = paste(prefix,"PCA plot of PC",x," vs PC",y,sep=""), las =1, 
       cex=plot.pt.cex, cex.axis=plot.cex-0.3, cex.lab=plot.cex-0.15, cex.main=plot.cex)
  abline(0,0, lty=5)
  abline(0,180, lty=5)
  for (sample in 2:number.indv) {
    points(pca1$li[sample,x], pca1$li[sample,y], col=pop.features[sample,1], pch=pop.features[sample,2], cex=plot.pt.cex)
  }
}

m=matrix(c(1,1,1,2,2,2,
           1,1,1,2,2,2,
           1,1,1,2,2,2,
           3,3,3,3,3,3), ncol =6, byrow = TRUE)

image.size=1
wid=600*image.size
hig=380*image.size
plot.size=2
point.size=2

png(file=paste("PCA_",outputfile.name,"_combined_font.png",sep=""), height = hig, width = wid, pointsize =12)
layout(m, heights = c(3,1))
par(mar=c(5, 5.5, 4, 2), mgp = c(3.5, 1, 0))
plot.PCA(1, 2, ld.pos.PC12, plot.cex=plot.size, plot.pt.cex=point.size-0.5, plot.ncol=1, prefix = "A) ")
plot.PCA(1, 3, ld.pos.PC12, plot.cex=plot.size, plot.pt.cex=point.size-0.5, plot.ncol=1, prefix = "B) ")
par(mar=c(1, 1, 1, 1))
plot.new()

legend(x="center", legend=c("P. × adulterinum", "P. gallicum", "P. nigrum", "P. pyrenaicum", "P. spicatum East", "P. spicatum West"), 
       col = unique(pop.features[[1]]), pch = unique(pop.features[[2]]), pt.cex=point.size, cex=plot.size, ncol=3, y.intersp=1, 
       bty='n', text.font = 1)
dev.off()
