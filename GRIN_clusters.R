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
excluded <- c(names(rowtotals[rowtotals < 10]), 'District_of_Columbia')
#preplotdata <- preplotdata/(coltotals)*100
#preplotdata <- as.data.frame(t(t(preplotdata)/(rowtotals)*100)) #totals by plot
preplotdata <- preplotdata[order(row.names(preplotdata), decreasing = FALSE),sort(colnames(preplotdata), decreasing = FALSE)]
preplotdata <- preplotdata[,!(names(preplotdata) %in% excluded)]
plotdata <- t(preplotdata)
#betasim
betasim <- function(p){
  d <- matrix(1, nrow = nrow(p), ncol = nrow(p))
  rownames(d) <- rownames(p)
  colnames(d) <- rownames(p)
  for(j in 1:nrow(p)){
    for(k in 1:nrow(p)){
      d[j,k] <- 1-sum((p[j,]*p[k,])^0.5)/min(sum(p[j,]), sum(p[k,]))
    }}
  d<-as.dist(d)
}
betasim2 <- function(p){
  d <- matrix(1, nrow = nrow(p), ncol = nrow(p))
  rownames(d) <- rownames(p)
  colnames(d) <- rownames(p)
  for(j in 1:nrow(p)){
    for(k in 1:nrow(p)){
      d[j,k] <- 1-sum((p[j,]*p[k,])^0.5)/sqrt(sum(p[j,])*sum(p[k,]))
    }}
  d<-as.dist(d)
}

#makeplot

makeplot <- function(a,d,t,k){
  filename <- paste0('output/GRIN_',a,'.png')
  t <- as.hclust(t)
  #make cuts and reformat dendrogram
  ngroups=k
  groups <- cutree(t, k = ngroups)
  
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
  
  newlabels <- t$labels
  newlabels <- as.data.frame(newlabels)
  newlabels$row <- row(newlabels)
  newlabels <- merge(newlabels, groupdf, by.x='newlabels', by.y ='soilplot')
  newlabels$newlabels <- paste(newlabels$clust, newlabels$newlabels)
  newlabels <- newlabels[order(newlabels$row),1]
  newtree <- t
  newtree$labels <- newlabels
  
  dend1 <- color_branches(as.hclust(newtree), k = ngroups)
  dend1 <- color_labels(dend1, k = ngroups)
  
  #output file
  
  w <- 800
  h <- nrow(plotdata)*12+80
  u <- 12
  png(filename=filename,width = w, height = h, units = "px", pointsize = u)
  
  par(mar = c(2,0,1,13))
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', a,'method of', 'GRIN genera'), font=1, cex=0.84)
  #rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
  dev.off()
  
}

#analysis method ###Important to make sure all dist matrices are as.dist so that cluster analysis interprets correctly.###
a <- 'bray-agnes' 
if (T){
  a1 <- 'bray-agnes'
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t1 <- agnes(d, method='average')
  makeplot(a1,d,t1,k)
}
if (T){
  a2 <- 'jaccard-agnes'
  k=14
  d <- ((vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)))
  t2 <- agnes(d, method='average')
  makeplot(a2,d,t2,k)
}
if (T){
  a3 <- 'bray-ward' 
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t3 <- agnes(d, method='ward')
  makeplot(a3,d,t3,k)
}
if (T){
  a4 <- 'bray-diana' 
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t4 <- diana(d)
  makeplot(a4,d,t4,k)
}
if (T){
  a5 <- 'simpson-agnes'
  k=14
  d <- as.dist(simil(plotdata,method='Simpson'))
  t5 <- agnes(d, method='average')
  makeplot(a5,d,t5,k)
}
if (T){
  a6 <- 'bray-weighted'
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t6 <- agnes(d, method='weighted')
  makeplot(a6,d,t6,k)
}
if (F){
  a <- 'bray-single'
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t <- agnes(d, method='single')
  makeplot(a,d,jactree,k)
}
if (F){
  a <- 'bray-complete' 
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t <- agnes(d, method='complete')
  makeplot(a,d,jactree,k)
}
if (F){
  a <- 'simpson-diana'
  k=14
  d <- ((simil(plotdata,method='Simpson')))
  t <- diana(d)
  makeplot(a,d,jactree,k)
}
if (F){
  a <- 'simpson-ward'
  k=14
  d <- ((simil(plotdata,method='Simpson')))
  t <- agnes(d, method='ward')
  makeplot(a,d,jactree,k)
}
if (F){
  a <- 'euclid-ward'
  k=14
  t <- agnes(plotdata, method='ward')
  makeplot(a,d,jactree,k)
}
if (F){
  a <- 'agnes-betasim'
  k=14
  d <- betasim2(plotdata)
  t <- agnes(distbet, method='average')
  makeplot(a,distbet,tree,k)
}
if (F){
  a <- 'ward-betasim'
  k=14
  d <- betasim2(plotdata)
  t <- agnes(distbet, method='ward')
  makeplot(a,distbet,tree,k)
}
if (F){
  a <- 'diana-betasim'
  k=14
  d <- betasim2(plotdata)
  t <- diana(distbet)
  makeplot(a,distbet,tree,k)
}

if (T){
  a1 <- 'agnes-kulczynski'
  k=14
  d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
  t1 <- agnes(d, method='average')
  makeplot(a1,d,t1,k)
}
if (T){
  a1 <- 'ward-kulczynski'
  k=14
  d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
  t1 <- agnes(d, method='ward')
  makeplot(a1,d,t1,k)
}


#constraint compare
distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
distbet <- betasim2(plotdata)
distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
constdist <- read.delim("data/GRIN/constraints.txt")
rownames(constdist) <- constdist[,1]
constdist <- constdist[,-1]
constdist <- constdist[order(row.names(constdist), decreasing = FALSE),sort(colnames(constdist), decreasing = FALSE)]
constdist<- constdist[!rownames(constdist) %in% excluded, !names(constdist) %in% excluded]


 
coph.agnes <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distbray, method='average'))))
coph.jaccard <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distjac, method='average'))))
coph.dianajaccard <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(diana(distjac))))
coph.single <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distbray, method='single'))))
coph.complete <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distbray, method='complete'))))
coph.ward <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distbray, method='ward'))))
coph.diana <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(diana(distbray))))
coph.sim <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distsim, method='average'))))
coph.wardbet <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distbet, method='ward'))))
coph.dianabet <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(diana(distbet))))
coph.agnesbet <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distbet, method='average'))))
coph.kulc <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distkulc, method='average'))))
coph.wardkulc <- cor(cophenetic(as.hclust(agnes(distsim, method='average'))), cophenetic(as.hclust(agnes(distkulc, method='ward'))))


coph.agnes
coph.jaccard
coph.dianajaccard
coph.ward
coph.diana
coph.single
coph.complete
coph.sim
coph.agnesbet
coph.wardbet
coph.dianabet
coph.kulc
coph.wardkulc
#validity using silhouette index
#----
distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)

tbrayagnes <- distbray %>% agnes(method = 'average')
tbrayjac <- distjac %>% agnes(method = 'average')
tbrayward <- distbray %>% agnes(method = 'ward')
tbraydiana <- distbray %>% diana 
tsimpagnes <- distsim %>% agnes(method = 'average') 
tkulcagnes <- distkulc %>% agnes(method = 'average')
tkulcward <- distkulc %>% agnes(method = 'ward')
k <- 2
klevel <- 0
sil.bray <- 0
sil.jac <- 0
sil.sim <- 0
sil.ward <- 0
sil.diana <- 0
sil.kmeans <- 0
sil.kulc <- 0
sil.wardkulc <- 0
for (k in 2:20){
sil.bray1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
sil.jac1 <- (tbrayjac %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
sil.sim1 <- (tsimpagnes %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
sil.ward1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
sil.diana1 <- (tbraydiana %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
sil.kmeans1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distkulc))[,3] %>% mean
sil.kulc1 <- (tkulcagnes %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
sil.wardkulc1 <- (tkulcward %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean


klevel <- c(klevel, k)
sil.bray <- c(sil.bray, sil.bray1)
sil.jac <- c(sil.jac, sil.jac1)
sil.sim <- c(sil.sim, sil.sim1)
sil.ward <- c(sil.ward, sil.ward1)
sil.diana <- c(sil.diana, sil.diana1)
sil.kmeans <- c(sil.kmeans, sil.kmeans1)
sil.kulc <- c(sil.kulc, sil.kulc1)
sil.wardkulc <- c(sil.wardkulc, sil.wardkulc1)
}
sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.kulc,sil.wardkulc))
sil.table <- sil.table[-1,]


#indicator analysis
k=15
distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
groups <- distbray %>% agnes(method = 'average') %>% cutree(k=k) 
tree <- agnes(distbray, method='average')
a <- 'test2'
#makeplot(a, distbray,tree,k)






spp.freq <- syntable(plotdata, groups)
spp.freq <- spp.freq$syntable
spp.mean <- syntable(plotdata, groups,  type = "mean")
spp.mean <- spp.mean$syntable
#reanalysis after grouping
k = 5
tspp.freq <- t(spp.freq)

if (T){
  a <- 'spp-freq-groups' 
  jdist <- as.data.frame(as.matrix(vegdist(tspp.freq, method='bray', binary=FALSE, na.rm=T)))
  tre <- agnes(jacdist, method='average')
  makeplot(a,jdist,tre,k)
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
#spp.diff <- syntable(plotdata, groups,  type = "diffspec")
#spp.diffonly <- rownames(spp.diff$onlydiff)
#spp.diff <- spp.diff$syntable
#spp.freqdif <- subset(spp.freq, rownames(spp.freq) %in% spp.diffonly)
#spp.diffonly <- subset(spp.diff, rownames(spp.diff) %in% spp.diffonly)




silanalysis2 <- function(input){
  distbray <- vegdist(input, method='bray', binary=FALSE, na.rm=T)
  distjac <- vegdist(input, method='jaccard', binary=FALSE, na.rm=T)
  distsim <- as.dist(simil(input,method='Simpson'))

  maxcluster <- min(20, nrow(input)-1)
  k <- 2
  klevel <- 0
  sil.bray <- 0
  sil.jac <- 0
  sil.sim <- 0
  sil.ward <- 0
  sil.diana <- 0
  
  for (k in 2:maxcluster){
    sil.bray1 <- (distbray %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    sil.bray1 <- (aggregate(sil.bray1[,2], by=list(sil.bray1[,1]), FUN='mean'))[,2] %>% mean()
    sil.jac1 <- (distjac %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    sil.jac1 <- (aggregate(sil.jac1[,2], by=list(sil.jac1[,1]), FUN='mean'))[,2] %>% mean()
    sil.sim1 <- (distsim %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    sil.sim1 <- (aggregate(sil.sim1[,2], by=list(sil.sim1[,1]), FUN='mean'))[,2] %>% mean()
    sil.ward1 <- (distbray %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbray))[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    sil.ward1 <- (aggregate(sil.ward1[,2], by=list(sil.ward1[,1]), FUN='mean'))[,2] %>% mean()
    sil.diana1 <- (distbray %>% diana %>% cutree(k=k) %>% silhouette(distbray))[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    sil.diana1 <- (aggregate(sil.diana1[,2], by=list(sil.diana1[,1]), FUN='mean'))[,2] %>% mean()
    
    klevel <- c(klevel, k)
    sil.bray <- c(sil.bray, sil.bray1)
    sil.jac <- c(sil.jac, sil.jac1)
    sil.sim <- c(sil.sim, sil.sim1)
    sil.ward <- c(sil.ward, sil.ward1)
    sil.diana <- c(sil.diana, sil.diana1)
  }
  sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana))
  sil.table <- sil.table[-1,]
  sil.table2 <<- sil.table
  return(sil.table)}
silanalysis2(plotdata)
##### hybrid analysis ----
### upgma to ward
k <- min(max(floor(nrow(plotdata)/10),2),10)
d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
t <- agnes(d, method='average')
groups <- cutree(t, k = k)
name <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(name, clust))
groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
maxcluster <- max(groupdf$clust)
numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
whichrecords <- which(groupdf$clust == 0)
if (nrow(groupdf[groupdf$clust == 0,]) != 0){
  for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
    groupdf[whichrecords[i],]$clust <- maxcluster+i}}
groupdf$p <- 1/2^0.5
groupdf$clust <- paste0('c',groupdf$clust)
m <- makecommunitydataset(groupdf, row = 'name', column = 'clust', value = 'p')
d1 <- vegdist(m, method = 'euclidean')
d1tab <- as.matrix(d1)
dtab <- as.matrix(d)
d2 <- d + d1/20
t1 <- agnes(d2, method='ward')
k2 <- 14
a1 <- 'hybrid_upgma_to_ward'
makeplot(a1,d,t1,k2)
### ward to upgma
k <- min(max(floor(nrow(plotdata)/10),50),nrow(plotdata))
d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
t <- agnes(d, method='ward')
groups <- cutree(t, k = k)
name <- names(groups)
clust <- unname(groups)
groupdf <- as.data.frame(cbind(name, clust))
groupdf$clust <- (as.numeric(as.character(groupdf$clust)))
maxcluster <- max(groupdf$clust)
numberzeros <- nrow(groupdf[(groupdf$clust == 0),])
whichrecords <- which(groupdf$clust == 0)
if (nrow(groupdf[groupdf$clust == 0,]) != 0){
  for (i in 1:numberzeros){ #assign all zero clusters to unique cluster number.
    groupdf[whichrecords[i],]$clust <- maxcluster+i}}
groupdf$p <- 1/2^0.5
groupdf$clust <- paste0('c',groupdf$clust)
m <- makecommunitydataset(groupdf, row = 'name', column = 'clust', value = 'p')
d1 <- vegdist(m, method = 'euclidean')
d1tab <- as.matrix(d1)
dtab <- as.matrix(d)
d2 <- d + d1/20
t1 <- agnes(d2, method='average')
k2 <- 14
a1 <- 'hybrid_ward_to_upgma'
makeplot(a1,d,t1,k2)
#### cophenetic bridge----
k2 <-14
d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
tw <- agnes(d, method='ward')
tu <- agnes(d, method='average')
makeplot('testupgma',d,tu,k2)
length(d)
makeplot('testward',d,tw,k2)

coph.tw <- cophenetic(tw)/mean(cophenetic(tw))
d2 <- as.matrix(coph.tw)
d2a <- as.dist(ifelse(d2 > 1,1,d2))
t2 <- agnes(d2a, method='average')
makeplot('test2',d,t2,k2)

coph.tu <- cophenetic(tu)/mean(cophenetic(tu))
d2 <- as.matrix(coph.tu)
d2b <- as.dist(ifelse(d2 < 1,1,d2))
t2 <- agnes(d2b, method='average')
makeplot('test2',d,t2,k2)

d3 <- (coph.tw^(1/2)+d^2*2)/3
t3 <- agnes(d3, method='average')
makeplot('testupgmaward',d,t3,k2)

#isopam----
library(isopam)
pamtree <- isopam(plotdata, distance = 'kulczynski')

if (T){
  a <- 'isopam-kulczynski' 
  k = 14
  d <- ((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
  t <- pamtree$dendro
  makeplot(a,d,t,k)
}
pamtree$analytics
table <- isotab(pamtree, level = 1)
table <- table$tab
as.numeric(as.character((table$`1`)))
table$x1 <- as.numeric(gsub("[^0-9.-]", "", (table$`1`)))
table$x2 <- as.numeric(gsub("[^0-9.-]", "", (table$`2`)))
table$dif <- table$x1-table$x2
coph <- cophenetic(t)
d2 <- (coph/mean(coph)*2 + d/mean(d))/3
t2 <- agnes(d2, method = 'average')
makeplot('kulczynski-isopam-agnes-hybrid',d2,t2,k)



distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)

tbrayagnes <- distbray %>% agnes(method = 'average')
tbrayjac <- distjac %>% agnes(method = 'average')
tbrayward <- distbray %>% agnes(method = 'ward')
tbraydiana <- distbray %>% diana 
tsimpagnes <- distsim %>% agnes(method = 'average') 
tkulcagnes <- distkulc %>% agnes(method = 'average')
tkulcward <- distkulc %>% agnes(method = 'ward')
k <- 2
klevel <- 0
sil.bray <- 0
sil.jac <- 0
sil.sim <- 0
sil.ward <- 0
sil.diana <- 0
sil.kmeans <- 0
sil.kulc <- 0
sil.wardkulc <- 0
sil.isopam <- 0
for (k in 2:20){
  sil.bray1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.jac1 <- (tbrayjac %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.sim1 <- (tsimpagnes %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.ward1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.diana1 <- (tbraydiana %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.kmeans1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distkulc))[,3] %>% mean
  sil.kulc1 <- (tkulcagnes %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.wardkulc1 <- (tkulcward %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.isopam1 <- (pamtree$dendro %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  
  klevel <- c(klevel, k)
  sil.bray <- c(sil.bray, sil.bray1)
  sil.jac <- c(sil.jac, sil.jac1)
  sil.sim <- c(sil.sim, sil.sim1)
  sil.ward <- c(sil.ward, sil.ward1)
  sil.diana <- c(sil.diana, sil.diana1)
  sil.kmeans <- c(sil.kmeans, sil.kmeans1)
  sil.kulc <- c(sil.kulc, sil.kulc1)
  sil.wardkulc <- c(sil.wardkulc, sil.wardkulc1)
  sil.isopam <- c(sil.isopam, sil.isopam1)
}
sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.kulc,sil.wardkulc,sil.isopam))
sil.table <- sil.table[-1,]


cor(cophenetic(tbrayagnes), cophenetic(pamtree$dendro))
cor(cophenetic(tbrayjac), cophenetic(pamtree$dendro))
cor(cophenetic(tbraydiana), cophenetic(pamtree$dendro))
cor(cophenetic(tbrayward), cophenetic(pamtree$dendro))
cor(cophenetic(tsimpagnes), cophenetic(pamtree$dendro))
cor(cophenetic(tkulcagnes), cophenetic(pamtree$dendro))
cor(cophenetic(tkulcward), cophenetic(pamtree$dendro))

