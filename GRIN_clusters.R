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
library(proxy)
library(goeveg)
#import
preplotdata <- read.delim("data/GRIN/altgeogrin.txt")
rownames(preplotdata) <- preplotdata[,1]
preplotdata <- preplotdata[,-1]
excluded <- ''
#excluded <- c('Alabama', 'Arizona', 'Arkansas', 'British_Columbia', 'Connecticut', 'Delaware', 'District_of_Columbia', 'Georgia', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kentucky', 'Louisiana', 'Maryland', 'Massachusetts', 'Minnesota', 'Mississippi', 'Montana', 'Nebraska', 'Nevada', 'New_Brunswick', 'New_Hampshire', 'New_Jersey', 'New_Mexico', 'New_York', 'Newfoundland', 'North_Carolina', 'Northwest_Territory', 'Nova_Scotia', 'Ohio', 'Oklahoma', 'Ontario', 'Oregon', 'Pennsylvania', 'Prince_Edward_Island', 'Quebec', 'Rhode_Island', 'Saskatchewan', 'South_Dakota', 'Vermont', 'Virginia', 'Wisconsin', 'Wyoming', 'Yukon_Territory')
coltotals <- (apply(preplotdata, MARGIN = 1, FUN = 'sum' ))
rowtotals <- (apply(preplotdata, MARGIN = 2, FUN = 'sum' ))
#excluded <- c(names(rowtotals[rowtotals < 200]), 'District_of_Columbia')
#preplotdata <- preplotdata/(coltotals)*100
#preplotdata <- as.data.frame(t(t(preplotdata)/(rowtotals)*100)) #totals by plot
preplotdata <- preplotdata[order(row.names(preplotdata), decreasing = FALSE),sort(colnames(preplotdata), decreasing = FALSE)]
preplotdata <- preplotdata[,!(names(preplotdata) %in% excluded)]
plotdata <- t(preplotdata)


makeplot <- function(amethod,jacdist,jactree,k){
  filename <- paste0('output/GRIN_',amethod,'.png')
  
  #make cuts and reformat dendrogram
  ngroups=k
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
  
}

#analysis method ###Important to make sure all dist matrices are as.dist so that cluster analysis interprets correctly.###
amethod <- 'bray-agnes' 
if (T){
  amethod <- 'bray-agnes'
  k=15
  jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='average')
  makeplot(amethod,jacdist,jactree,k)
}
if (T){
  amethod <- 'jaccard-agnes'
  k=15
  jacdist <- ((vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='average')
  makeplot(amethod,jacdist,jactree,k)
}
if (T){
  amethod <- 'bray-ward' 
  k=9
  jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='ward')
  makeplot(amethod,jacdist,jactree,k)
}
if (T){
  amethod <- 'bray-diana' 
  k=10
  jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- diana(jacdist)
  makeplot(amethod,jacdist,jactree,k)
}
if (T){
  amethod <- 'bray-single'
  k=15
  jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='single')
  makeplot(amethod,jacdist,jactree,k)
}
if (T){
  amethod <- 'bray-complete' 
  k=7
  jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(jacdist, method='complete')
  makeplot(amethod,jacdist,jactree,k)
}

if (T){
  amethod <- 'simpson-agnes'
  k=4
  jacdist <- as.dist(simil(plotdata,method='Simpson'))
  jactree <- agnes(jacdist, method='average')
  makeplot(amethod,jacdist,jactree,k)
}
if (F){
  amethod <- 'simpson-diana'
  k=4
  jacdist <- ((simil(plotdata,method='Simpson')))
  jactree <- diana(jacdist)
  makeplot(amethod,jacdist,jactree,k)
}
if (F){
  amethod <- 'simpson-ward'
  k=4
  jacdist <- ((simil(plotdata,method='Simpson')))
  jactree <- agnes(jacdist, method='ward')
  makeplot(amethod,jacdist,jactree,k)
}
if (F){
  amethod <- 'euclid-ward'
  k=5
  jactree <- agnes(plotdata, method='ward')
  makeplot(amethod,jacdist,jactree,k)
}


library(phytools)

if (F){
  amethod <- 'bray-nj' 
  jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- nj(jacdist)
  filename <- paste0('output/GRIN_',amethod,'.png')
  jactree <- root(jactree, outgroup = 'North.Central_Pacific')
  w <- 800
  h <- nrow(plotdata)*12+80
  u <- 12
  
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  par(mar = c(2,0,1,13))
  plot(((jactree)), main=paste('floristic simularity', amethod,'method of', 'GRIN genera'), label.offset=0.05, direction='right', font=1, cex=0.85)
  dev.off()
  filename <- paste0('output/GRIN_',amethod,'-ultrametric.png')
  jactree3 <- force.ultrametric(jactree, method=c("extend"))
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  par(mar = c(2,0,1,13))
  plot(((jactree3)), main=paste('floristic simularity', amethod,'method of', 'GRIN genera-ultrametric'), label.offset=0.05, direction='right', font=1, cex=0.85)
  dev.off()
  
  
  write.csv(jacdist, 'output/braydist.csv')
}

#constraint compare
distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
constdist <- read.delim("data/GRIN/constraints.txt")
rownames(constdist) <- constdist[,1]
constdist <- constdist[,-1]
constdist <- constdist[order(row.names(constdist), decreasing = FALSE),sort(colnames(constdist), decreasing = FALSE)]
constdist<- constdist[!rownames(constdist) %in% excluded, !names(constdist) %in% excluded]


 
coph.agnes <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='average'))))
coph.jaccard <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distjac, method='average'))))
coph.dianajaccard <- cor(as.dist(constdist), cophenetic(as.hclust(diana(distjac))))
coph.single <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='single'))))
coph.complete <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='complete'))))
coph.ward <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='ward'))))
coph.diana <- cor(as.dist(constdist), cophenetic(as.hclust(diana(distbray))))
coph.sim <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distsim, method='average'))))
coph.wardsim <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distsim, method='ward'))))
coph.dianasim <- cor(as.dist(constdist), cophenetic(as.hclust(diana(distsim))))


coph.agnes
coph.jaccard
coph.dianajaccard
coph.ward
coph.diana
coph.single
coph.complete
coph.sim
coph.wardsim
coph.dianasim
#validity using silhouette index
#----
distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
tree <- agnes(distbray, method='average')
cuttree <- cutree(tree, k = 8)
siltree <- silhouette(cuttree, distbray)

summary(siltree)
mean(siltree[,3])
k <- 2
klevel <- 0
sil.bray <- 0
sil.jac <- 0
sil.sim <- 0
sil.ward <- 0
sil.diana <- 0
sil.kmeans <- 0
sil.single <- 0
sil.complete <- 0
sil.wardeuc <- 0
sil.kmeanseuc <- 0
for (k in 2:20){
sil.bray1 <- (distbray %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.jac1 <- (distjac %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.sim1 <- (distsim %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.ward1 <- (distbray %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.diana1 <- (distbray %>% diana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.kmeans1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
sil.single1 <- (distbray %>% agnes(method = 'single') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.complete1 <- (distbray %>% agnes(method = 'complete') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.wardeuc1 <- (plotdata %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
sil.kmeanseuc1 <- (kmeans(plotdata, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean

klevel <- c(klevel, k)
sil.bray <- c(sil.bray, sil.bray1)
sil.jac <- c(sil.jac, sil.jac1)
sil.sim <- c(sil.sim, sil.sim1)
sil.ward <- c(sil.ward, sil.ward1)
sil.diana <- c(sil.diana, sil.diana1)
sil.kmeans <- c(sil.kmeans, sil.kmeans1)
sil.single <- c(sil.single, sil.single1)
sil.complete <- c(sil.complete, sil.complete1)
sil.wardeuc <- c(sil.wardeuc, sil.wardeuc1)
sil.kmeanseuc <- c(sil.kmeanseuc, sil.kmeanseuc1)
}
sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.single,sil.complete))
sil.table <- sil.table[-1,]


#indicator analysis
k=15
distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
groups <- distbray %>% agnes(method = 'average') %>% cutree(k=k) 
tree <- agnes(distbray, method='average')
amethod <- 'test2'
#makeplot(amethod, distbray,tree,k)






spp.freq <- syntable(plotdata, groups)
spp.freq <- spp.freq$syntable
spp.mean <- syntable(plotdata, groups,  type = "mean")
spp.mean <- spp.mean$syntable
#reanalysis after grouping
k = 5
tspp.freq <- t(spp.freq)

if (T){
  amethod <- 'spp-freq-groups' 
  jdist <- as.data.frame(as.matrix(vegdist(tspp.freq, method='bray', binary=FALSE, na.rm=T)))
  tre <- agnes(jacdist, method='average')
  makeplot(amethod,jdist,tre,k)
}



#----
distbray2 <- as.data.frame(as.matrix(distbray))
distbray2 <- (distbray2+constdist)/2
makeplot('test',distbray,agnes(distbray, method='average'))
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
