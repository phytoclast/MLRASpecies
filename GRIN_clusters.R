library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
#library(plyr)
library(dplyr)
library(dynamicTreeCut)
library(rpart)
library(rpart.plot)
library(goeveg)
#import
plotdata <- read.delim("data/GRIN/altgeogrin.txt")
rownames(plotdata) <- plotdata[,1]
plotdata <- plotdata[,-1]
plotdata <- t(plotdata)
#analysis method
amethod <- 'bray-agnes' 
if (T){
  amethod <- 'bray-agnes' 
  jacdist <- as.data.frame(as.matrix(vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='average')
}
if (F){
  amethod <- 'bray-single' 
  jacdist <- as.data.frame(as.matrix(vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='single')
}
if (F){
  amethod <- 'bray-complete' 
  jacdist <- as.data.frame(as.matrix(vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='complete')
}
if (F){
  amethod <- 'bray-diana' 
  jacdist <- as.data.frame(as.matrix(vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- diana(jacdist)
}
if (F){
  amethod <- 'bray-ward' 
  jacdist <- as.data.frame(as.matrix(vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='ward')
}
if (T){
  amethod <- 'jaccard-agnes' 
  jacdist <- as.data.frame(as.matrix(vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='average')
}
filename <- paste0('output/GRIN_',amethod,'.png')

#make cuts and reformat dendrogram
ngroups=8
groups <- cutree(jactree, k = ngroups)

soilplot <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(soilplot, clust))
groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
maxcluster <- max(groupdf$clust)
numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
whichrecords <- which(groupdf$clust == 0)
if (nrow(groupdf[groupdf$clust == 0,]) != 0){
  for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
    groupdf[whichrecords[i],]$clust <- maxcluster+i}}

newlabels <- jactree$order.lab
newlabels <- as.data.frame(newlabels)
newlabels$row <- row(newlabels)
newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
newlabels <- newlabels[order(newlabels$row),1]
newtree <- jactree
newtree$order.lab <- newlabels

dend1 <- color_branches(as.hclust(newtree), k = ngroups)
dend1 <- color_labels(dend1, k = ngroups)

#output file

w <- 800
h <- nrow(plotdata)*12+80
u <- 12
png(filename=filename,width = w, height = h, units = "px", pointsize = u)

par(mar = c(2,0,1,13))
plot(dend1, horiz = TRUE, main=paste('floristic simularity', amethod,'method of', 'GRIN genera'), font=1, cex=0.84)
#rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
dev.off()

#indicator species analysis
#----
spp.freq <- syntable(plotdata, groups)
spp.freq <- spp.freq$syntable
spp.mean <- syntable(plotdata, groups,  type = "mean")
spp.mean <- spp.mean$syntable
spp.med <- syntable(plotdata, groups,  type = "median")
spp.med <- spp.med$syntable
spp.diff <- syntable(plotdata, groups,  type = "diffspec")
spp.diffonly <- rownames(spp.diff$onlydiff)
spp.diff <- spp.diff$syntable
spp.freqdif <- subset(spp.freq, rownames(spp.freq) %in% spp.diffonly)
spp.diffonly <- subset(spp.diff, rownames(spp.diff) %in% spp.diffonly)
