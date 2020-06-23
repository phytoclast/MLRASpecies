# load required libraries
library(BiodiversityR)
library(cluster)
library(ape)
library(dendextend)
library(dplyr)
library(dynamicTreeCut)
library(rpart)
library(rpart.plot)
library(goeveg)
library(proxy)
library(goeveg)
library(foreign)
######################################----

#subsections
plotdata <- readRDS('data/usfssubsecmatrix.RDS')
distbray <- readRDS('data/usfssubsec.distbray.RDS')
distjac <- readRDS('data/usfssubsec.distjac.RDS')
distsim <- readRDS('data/usfssubsec.distsim.RDS')
distkulc <- readRDS('data/usfssubsec.distkulc.RDS')

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


#validity using silhouette index
#----
#distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
#distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
#distsim <- as.dist(simil(plotdata,method='Simpson'))
#distbeta <- betasim(plotdata)
#distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
#distbeta2 <- betasim2(plotdata)
#saveRDS(distbray,'data/usfssubsec.distbray.RDS')
#saveRDS(distjac,'data/usfssubsec.distjac.RDS')
#saveRDS(distsim,'data/usfssubsec.distsim.RDS')
#saveRDS(distbeta,'data/usfssubsec.distbeta.RDS')
#saveRDS(distkulc,'data/usfssubsec.distkulc.RDS')
#saveRDS(distbeta2,'data/usfssubsec.distbeta2.RDS')

#### test different algorythms ----
t.bray <- distbray %>% agnes(method = 'average')
t.jac <- distjac %>% agnes(method = 'average') 
t.ward <- distbray %>% agnes(method = 'ward') 
t.diana <- distbray %>% diana 
t.kulc <- distkulc %>% agnes(method = 'average') 
t.wkul <- distkulc %>% agnes(method = 'ward')
d.ku <- cophenetic(t.kulc)/mean(cophenetic(t.kulc))
d.kw <- cophenetic(t.wkul)/mean(cophenetic(t.wkul))
d.hybrid <- (d.kw^0.5+distkulc^2)/2
t.hybrid <- d.hybrid %>% agnes(method = 'average')

k <- 2
klevel <- 0
sil.bray <- 0
sil.jac <- 0
sil.beta <- 0
sil.ward <- 0
sil.diana <- 0
sil.kmeans <- 0
sil.jdiana <- 0
sil.complete <- 0
sil.wardeuc <- 0
sil.kmeanseuc <- 0
sil.kulc <- 0
sil.wardkulc <- 0
sil.hybrid <- 0



for (k in 2:24){
  sil.bray1 <- (t.bray %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.jac1 <- (t.jac %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.ward1 <- (t.ward %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.diana1 <- (t.diana %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.kulc1 <- (t.kulc %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.wardkulc1 <- (t.wkul %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  sil.hybrid1 <- (t.hybrid %>% cutree(k=k) %>% silhouette(distkulc))[,3]%>% mean
  

  klevel <- c(klevel, k)
  sil.bray <- c(sil.bray, sil.bray1)
  sil.jac <- c(sil.jac, sil.jac1)
  sil.ward <- c(sil.ward, sil.ward1)
  sil.diana <- c(sil.diana, sil.diana1)
  sil.kulc <- c(sil.kulc, sil.kulc1)
  sil.wardkulc <- c(sil.wardkulc, sil.wardkulc1)
  sil.hybrid <- c(sil.hybrid, sil.hybrid1)
}
sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.ward,sil.diana,sil.kulc,sil.wardkulc,sil.hybrid))
sil.table <- sil.table[-1,]
#saveRDS(sil.table,'output/usfssubsec.sil.table.RDS')


k.agnes <- c(2,3,4,8,10,16,24)
k.ward <- c(2,3,4,6,11,16,24)
k.diana <- c(2,3,4,7,10,16,24)



agnes.cut1<- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[1])
ward.cut1 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[1])
diana.cut1 <- distbray %>% diana %>% cutree(k=k.diana[1])
agnes.cut2 <- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[2])
ward.cut2 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[2])
diana.cut2 <- distbray %>% diana %>% cutree(k=k.diana[2])
agnes.cut3<- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[3])
ward.cut3 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[3])
diana.cut3 <- distbray %>% diana %>% cutree(k=k.diana[3])
agnes.cut4<- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[4])
ward.cut4 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[4])
diana.cut4 <- distbray %>% diana %>% cutree(k=k.diana[4])
agnes.cut5<- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[5])
ward.cut5 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[5])
diana.cut5 <- distbray %>% diana %>% cutree(k=k.diana[5])
agnes.cut6<- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[6])
ward.cut6 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[6])
diana.cut6 <- distbray %>% diana %>% cutree(k=k.diana[6])
agnes.cut7<- distbray %>% agnes(method = 'average') %>% cutree(k=k.agnes[7])
ward.cut7 <- distbray %>% agnes(method = 'ward') %>% cutree(k=k.ward[7])
diana.cut7 <- distbray %>% diana %>% cutree(k=k.diana[7])
simp.cut1<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[1])
simp.cut2<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[2])
simp.cut3<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[3])
simp.cut4<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[4])
simp.cut5<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[5])
simp.cut6<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[6])
simp.cut7<- distsim %>% agnes(method = 'average') %>% cutree(k=k.agnes[7])


usfsclust <- as.data.frame(cbind(agnes.cut1,agnes.cut2,agnes.cut3,agnes.cut4,agnes.cut5,agnes.cut6,agnes.cut7,
                                 ward.cut1,ward.cut2,ward.cut3,ward.cut4,ward.cut5,ward.cut6,ward.cut7,
                                 diana.cut1,diana.cut2,diana.cut3,diana.cut4,diana.cut5,diana.cut6,diana.cut7,
                                 simp.cut1,simp.cut2,simp.cut3,simp.cut4,simp.cut5,simp.cut6,simp.cut7))
usfsclust$subsect <- rownames(usfsclust)

kulc.cut1<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[1])
kulc.cut2<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[2])
kulc.cut3<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[3])
kulc.cut4<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[4])
kulc.cut5<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[5])
kulc.cut6<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[6])
kulc.cut7<- distkulc %>% agnes(method = 'average') %>% cutree(k=k.agnes[7])

kward.cut1<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[1])
kward.cut2<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[2])
kward.cut3<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[3])
kward.cut4<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[4])
kward.cut5<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[5])
kward.cut6<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[6])
kward.cut7<- distkulc %>% agnes(method = 'ward') %>% cutree(k=k.agnes[7])




usfsclust <- as.data.frame(cbind(kulc.cut1, kulc.cut1, kulc.cut2, kulc.cut3, kulc.cut4, kulc.cut5, kulc.cut6, kulc.cut7, kward.cut1, kward.cut2, kward.cut3, kward.cut4, kward.cut5, kward.cut6, kward.cut7))
usfsclust$subsect <- rownames(usfsclust)

write.dbf(usfsclust, 'output/usfsclust.dbf')

makeplot <- function(amethod,jacdist,jactree,k){
  filename <- paste0('output/USFS_',amethod,'.png')
  
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
  plot(dend1, horiz = TRUE, main=paste('floristic simularity', amethod,'method of', 'USFS subsections/BONAP Species'), font=1, cex=0.84)
  #rect.dendrogram(dend1, k = ngroups, horiz = TRUE)
  dev.off()
  
}

#analysis method ###Important to make sure all dist matrices are as.dist so that cluster analysis interprets correctly.###

if (T){
  amethod <- 'bray-agnes'
  k=10
  #jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- agnes(distbray, method='average')
  makeplot(amethod,distbray,tree,k)
}
if (T){
  amethod <- 'bray-ward' 
  k=11
  #jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  tree <- agnes(distbray, method='ward')
  makeplot(amethod,distbray,tree,k)
}
if (T){
  amethod <- 'bray-diana' 
  k=10
  #jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- diana(distbray)
  makeplot(amethod,distbray,tree,k)
}
if (T){
  amethod <- 'jaccard-diana' 
  k=10
  #jacdist <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  jactree <- diana(distjac)
  makeplot(amethod,distjac,tree,k)
}

#write.csv(as.matrix(distbray),'output/usfsdistbray.csv')
distsect <- read.delim('data/distsect.txt')
distprov <- read.delim('data/distprov.txt')
row.names(distprov) <- distprov[,1]
distprov <- as.dist(distprov[,-1])
row.names(distsect) <- distsect[,1]
distsect <- as.dist(distsect[,-1])
distboth <- (distprov+distsect)/2

constdist <- distprov

coph.agnes <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='average'))))
coph.jaccard <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distjac, method='average'))))
coph.dianajaccard <- cor(as.dist(constdist), cophenetic(as.hclust(diana(distjac))))
coph.single <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='single'))))
coph.complete <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='complete'))))
coph.ward <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distbray, method='ward'))))
coph.wardjaccard <- cor(as.dist(constdist), cophenetic(as.hclust(agnes(distjac, method='ward'))))
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
coph.wardjaccard



silanalysis2 <- function(input){
  #distbray <- vegdist(input, method='bray', binary=FALSE, na.rm=T)
  #distjac <- vegdist(input, method='jaccard', binary=FALSE, na.rm=T)
  #distsim <- as.dist(simil(input,method='Simpson'))
  
  maxcluster <- min(24, nrow(input)-1)
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
#grouping method
k <- 100#min(max(floor(nrow(plotdata)/10),2),10)
d <- distkulc #((vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)))
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
#merge method
t.kulc <- distkulc %>% agnes(method = 'average') 
t.wkul <- distkulc %>% agnes(method = 'ward')
d.ku <- cophenetic(t.kulc)/mean(cophenetic(t.kulc))
d.kw <- cophenetic(t.wkul)/mean(cophenetic(t.wkul))
d.hybrid <- (d.kw^0.5+distkulc^2)/2
t1 <- d.hybrid %>% agnes(method = 'average')



k.agnes <- c(2,3,4,8,10,16,24)
kulc.cut1<- t1 %>% cutree(k=k.agnes[1])
kulc.cut2<- t1 %>% cutree(k=k.agnes[2])
kulc.cut3<- t1 %>% cutree(k=k.agnes[3])
kulc.cut4<- t1 %>% cutree(k=k.agnes[4])
kulc.cut5<- t1 %>% cutree(k=k.agnes[5])
kulc.cut6<- t1 %>% cutree(k=k.agnes[6])
kulc.cut7<- t1 %>% cutree(k=k.agnes[7])

usfsclust <- as.data.frame(cbind(kulc.cut1, kulc.cut1, kulc.cut2, kulc.cut3, kulc.cut4, kulc.cut5, kulc.cut6, kulc.cut7))
usfsclust$subsect <- rownames(usfsclust)
write.dbf(usfsclust, 'output/usfsclust.dbf')
