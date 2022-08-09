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
library(optpart)
library(labdsv)
library(indicspecies)
#import
preplotdata <- read.delim("data/GRIN/altgeogrin.txt")
rownames(preplotdata) <- preplotdata[,1]
preplotdata <- preplotdata[,-1]
excluded <- ''
#excluded <- c('Alabama', 'Arizona', 'Arkansas', 'British_Columbia', 'Connecticut', 'Delaware', 'District_of_Columbia', 'Georgia', 'Idaho', 'Illinois', 'Indiana', 'Iowa', 'Kentucky', 'Louisiana', 'Maryland', 'Massachusetts', 'Minnesota', 'Mississippi', 'Montana', 'Nebraska', 'Nevada', 'New_Brunswick', 'New_Hampshire', 'New_Jersey', 'New_Mexico', 'New_York', 'Newfoundland', 'North_Carolina', 'Northwest_Territory', 'Nova_Scotia', 'Ohio', 'Oklahoma', 'Ontario', 'Oregon', 'Pennsylvania', 'Prince_Edward_Island', 'Quebec', 'Rhode_Island', 'Saskatchewan', 'South_Dakota', 'Vermont', 'Virginia', 'Wisconsin', 'Wyoming', 'Yukon_Territory')
coltotals <- (apply(preplotdata, MARGIN = 1, FUN = 'sum' ))
rowtotals <- (apply(preplotdata, MARGIN = 2, FUN = 'sum' ))
excluded <- c(names(rowtotals[rowtotals < 40]), 'District_of_Columbia')
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
a <- 'bray-flex' 
if (T){
  a1 <- 'bray-flex'
  k=14
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t1 <- flexbeta(d, beta= -0.25)
  makeplot(a1,d,t1,k)
}
a <- 'bray-flex3' 
if (T){
  a1 <- 'bray-flex25'
  k=8
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t1 <- flexbeta(d, beta= -0.25)
  makeplot(a1,d,t1,k)
}
a <- 'bray-flex1' 
if (T){
  a1 <- 'bray-flex15'
  k=3
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t1 <- flexbeta(d, beta= -0.15)
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
  k=4
  d <- ((vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)))
  t3 <- agnes(d, method='ward')
  makeplot(a3,d,t3,k)
}
if (T){
  a4 <- 'bray-diana' 
  k=2
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


distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata,method='Simpson'))
distkulc <- vegdist(plotdata, method='kulczynski', binary=FALSE, na.rm=T)
tbrayagnes <- distbray %>% agnes(method = 'average')
tbrayflex <- distbray %>% flexbeta(beta= -0.25)
tbrayflex1 <- distbray %>% flexbeta(beta= -0.1)
tbrayflex3 <- distbray %>% flexbeta(beta= -0.3)
tjacagnes <- distjac %>% agnes(method = 'average')
tbrayward <- distbray %>% agnes(method = 'ward')
tbraydiana <- distbray %>% diana 
tsimpagnes <- distsim %>% agnes(method = 'average') 
tkulcagnes <- distkulc %>% agnes(method = 'average')
tkulcward <- distkulc %>% agnes(method = 'ward')
k <- 2
klevel <- 0
sil.upgma <- 0
sil.flex <- 0
sil.flex1 <- 0
sil.flex3 <- 0
sil.opt <- 0
sil.jac <- 0
sil.sim <- 0
sil.ward <- 0
sil.diana <- 0
sil.kmeans <- 0
sil.kulc <- 0
sil.wardkulc <- 0

for (k in 2:20){#k=8
  sil.upgma.1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex.1 <- (tbrayflex %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex1.1 <- (tbrayflex1 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex3.1 <- (tbrayflex3 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  initialgrps <- tbrayward %>% cutree(k=k) 
  sil.opt.1 <- (optpart::optsil(initialgrps, distbray) %>% silhouette(distbray))[,3]%>% mean
  sil.jac.1 <- (tjacagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.sim.1 <- (tsimpagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.ward.1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.diana.1 <- (tbraydiana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.kmeans.1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
  sil.kulc.1 <- (tkulcagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.wardkulc.1 <- (tkulcward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  
  klevel <- c(klevel, k)
  sil.upgma <- c(sil.upgma, sil.upgma.1)
  sil.flex <- c(sil.flex, sil.flex.1)
  sil.flex1 <- c(sil.flex1, sil.flex1.1)
  sil.flex3 <- c(sil.flex3, sil.flex3.1)
  sil.opt <- c(sil.opt, sil.opt.1)
  sil.jac <- c(sil.jac, sil.jac.1)
  sil.sim <- c(sil.sim, sil.sim.1)
  sil.ward <- c(sil.ward, sil.ward.1)
  sil.diana <- c(sil.diana, sil.diana.1)
  sil.kmeans <- c(sil.kmeans, sil.kmeans.1)
  sil.kulc <- c(sil.kulc, sil.kulc.1)
  sil.wardkulc <- c(sil.wardkulc, sil.wardkulc.1)
  }
sil.table <- as.data.frame(cbind(klevel,sil.upgma,sil.flex,sil.flex1,sil.flex3,sil.opt,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.kulc,sil.wardkulc))
sil.table <- sil.table[-1,]
lis.table <- sil.table[,-1] %>% t() %>% as.data.frame()
lis.table <- lis.table %>% mutate(s1to8 = apply(lis.table[,1:8], MARGIN=1, FUN = 'mean'))
lis.table <- lis.table %>% mutate(s9to16 = apply(lis.table[,9:16], MARGIN=1, FUN = 'mean'))


cor(cophenetic(tbrayagnes), distbray)
cor(cophenetic(tbrayflex), distbray)
cor(cophenetic(tbrayflex1), distbray)
cor(cophenetic(tbrayflex3), distbray)
cor(cophenetic(tjacagnes), distbray)
cor(cophenetic(tbraydiana), distbray)
cor(cophenetic(tbrayward), distbray)
cor(cophenetic(tsimpagnes), distbray)
cor(cophenetic(tkulcagnes), distbray)
cor(cophenetic(tkulcward), distbray)

cor(cophenetic(tbrayagnes), cophenetic(tbrayward))
cor(cophenetic(tbrayflex), cophenetic(tbrayward))
cor(cophenetic(tbrayflex1), cophenetic(tbrayward))
cor(cophenetic(tbrayflex3), cophenetic(tbrayward))
cor(cophenetic(tjacagnes), cophenetic(tbrayward))
cor(cophenetic(tbraydiana), cophenetic(tbrayward))
cor(cophenetic(tbrayward), cophenetic(tbrayward))
cor(cophenetic(tsimpagnes), cophenetic(tbrayward))
cor(cophenetic(tkulcagnes), cophenetic(tbrayward))
cor(cophenetic(tkulcward), cophenetic(tbrayward))

#indval ----
plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]
initialgrps <- tbrayward %>% cutree(k=8) 

mp <- multipatt(plotdata1, initialgrps)
inds <- indicspecies::indicators(plotdata1, initialgrps)
cv <- coverage(plotdata1, mp)

oindval <- optpart::optimclass(plotdata1, initialgrps)
oindval$sums

#indval ----
plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]
indval.1 <- indval(plotdata1, initialgrps)
look <- as.data.frame(indval.1$indval)

totalinval <- mean(apply(look, MARGIN=1, FUN = 'mean'))

k <- 2
klevel <- 0
ind.upgma <- 0
ind.flex <- 0
ind.flex1 <- 0
ind.flex3 <- 0
ind.opt <- 0
ind.jac <- 0
ind.sim <- 0
ind.ward <- 0
ind.diana <- 0
ind.kmeans <- 0
ind.kulc <- 0
ind.wardkulc <- 0
for (k in 2:20){
  ind.upgma.1 <- tbrayagnes %>% cutree(k=k) 
  ind.flex.1 <- tbrayflex %>% cutree(k=k)
  ind.flex1.1 <- tbrayflex1 %>% cutree(k=k) 
  ind.flex3.1 <- tbrayflex3 %>% cutree(k=k) 
  initialgrps <- tbrayward %>% cutree(k=k) 
  ind.opt.1 <- optpart::optsil(initialgrps, distbray)
  ind.jac.1 <- tbrayjac %>% cutree(k=k)
  ind.sim.1 <- tsimpagnes %>% cutree(k=k)
  ind.ward.1 <- tbrayward %>% cutree(k=k)
  ind.diana.1 <- tbraydiana %>% cutree(k=k)
  ind.kmeans.1 <- kmeans(distbray, centers = k)$cluster 
  ind.kulc.1 <- tkulcagnes %>% cutree(k=k)
  ind.wardkulc.1 <- tkulcward %>% cutree(k=k) 
  ind.wardbraykulc.1 <- tkulcbrayward %>% cutree(k=k) 
  
  ind.upgma.1 <- indval(plotdata1, ind.upgma.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.flex.1 <- indval(plotdata1, ind.flex.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.flex1.1 <- indval(plotdata1, ind.flex1.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.flex3.1 <- indval(plotdata1, ind.flex3.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.opt.1 <- indval(plotdata1, ind.opt.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.jac.1 <- indval(plotdata1, ind.jac.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.sim.1 <- indval(plotdata1, ind.sim.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.ward.1 <- indval(plotdata1, ind.ward.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.diana.1 <- indval(plotdata1, ind.diana.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.kmeans.1 <- indval(plotdata1, ind.kmeans.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.kulc.1 <- indval(plotdata1, ind.kulc.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  ind.wardkulc.1 <- indval(plotdata1, ind.wardkulc.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  
  klevel <- c(klevel, k)
  ind.upgma <- c(ind.upgma, ind.upgma.1)
  ind.flex <- c(ind.flex, ind.flex.1)
  ind.flex1 <- c(ind.flex1, ind.flex1.1)
  ind.flex3 <- c(ind.flex3, ind.flex3.1)
  ind.opt <- c(ind.opt, ind.opt.1)
  ind.jac <- c(ind.jac, ind.jac.1)
  ind.sim <- c(ind.sim, ind.sim.1)
  ind.ward <- c(ind.ward, ind.ward.1)
  ind.diana <- c(ind.diana, ind.diana.1)
  ind.kmeans <- c(ind.kmeans, ind.kmeans.1)
  ind.kulc <- c(ind.kulc, ind.kulc.1)
  ind.wardkulc <- c(ind.wardkulc, ind.wardkulc.1)
}
ind.table <- as.data.frame(cbind(klevel,ind.upgma,ind.flex,ind.flex1,ind.flex3,ind.opt,ind.jac,ind.sim,ind.ward,ind.diana,ind.kmeans,ind.kulc,ind.wardkulc))
ind.table <- ind.table[-1,]
ind.table <- ind.table %>% mutate()
dni.table <- ind.table[,-1] %>% t() %>% as.data.frame()
dni.table <- dni.table %>% mutate(s1to8 = apply(dni.table[,1:8], MARGIN=1, FUN = 'mean'))
dni.table <- dni.table %>% mutate(s9to16 = apply(dni.table[,9:16], MARGIN=1, FUN = 'mean'))
write.csv(ind.table, 'output/ind.table.csv', row.names = F)


mean(apply(plotdata1, MARGIN = 2, FUN = 'sd'))
initialgrps <- tbrayward %>% cutree(k=8) 
groups <- sample(c(1:8), nrow(plotdata1), replace = T)

grptable <- aggregate(plotdata1, by=list(groups), FUN='mean') %>% .[,-1]
mean(apply(grptable, MARGIN = 2, FUN = 'sd'))
grptable <- aggregate(plotdata1, by=list(initialgrps), FUN='mean') %>% .[,-1]
mean(apply(grptable, MARGIN = 2, FUN = 'sd'))

#indval ----
plotdata.total <- apply(plotdata, MARGIN = 2, FUN = 'sum')
plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]

k <- 2
klevel <- 0
ind.upgma <- 0
ind.flex <- 0
ind.flex1 <- 0
ind.flex3 <- 0
ind.opt <- 0
ind.jac <- 0
ind.sim <- 0
ind.ward <- 0
ind.diana <- 0
ind.kmeans <- 0
ind.kulc <- 0
ind.wardkulc <- 0
for (k in 2:20){
    
  ind.upgma.0 <- tbrayagnes %>% cutree(k=k) 
  ind.flex.0 <- tbrayflex %>% cutree(k=k)
  ind.flex1.0 <- tbrayflex1 %>% cutree(k=k) 
  ind.flex3.0 <- tbrayflex3 %>% cutree(k=k) 
  initialgrps <- tbrayward %>% cutree(k=k) 
  ind.opt.0 <- optpart::optsil(initialgrps, distbray)
  ind.jac.0 <- tbrayjac %>% cutree(k=k)
  ind.sim.0 <- tsimpagnes %>% cutree(k=k)
  ind.ward.0 <- tbrayward %>% cutree(k=k)
  ind.diana.0 <- tbraydiana %>% cutree(k=k)
  ind.kmeans.0 <- kmeans(distbray, centers = k)$cluster 
  ind.kulc.0 <- tkulcagnes %>% cutree(k=k)
  ind.wardkulc.0 <- tkulcward %>% cutree(k=k) 
  
  ind.upgma.1 <- aggregate(plotdata1, by=list(ind.upgma.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.flex.1 <- aggregate(plotdata1, by=list(ind.flex.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.flex1.1 <- aggregate(plotdata1, by=list(ind.flex1.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.flex3.1 <- aggregate(plotdata1, by=list(ind.flex3.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.opt.1 <- aggregate(plotdata1, by=list(ind.opt.0$clustering), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.jac.1 <- aggregate(plotdata1, by=list(ind.jac.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.sim.1 <- aggregate(plotdata1, by=list(ind.sim.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.ward.1 <- aggregate(plotdata1, by=list(ind.ward.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.diana.1 <- aggregate(plotdata1, by=list(ind.diana.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.kmeans.1 <- aggregate(plotdata1, by=list(ind.kmeans.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.kulc.1 <- aggregate(plotdata1, by=list(ind.kulc.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  ind.wardkulc.1 <- aggregate(plotdata1, by=list(ind.wardkulc.0), FUN='mean') %>% .[,-1] %>% apply(MARGIN = 2, FUN = 'sd') %>% mean()
  
  klevel <- c(klevel, k)
  ind.upgma <- c(ind.upgma, ind.upgma.1)
  ind.flex <- c(ind.flex, ind.flex.1)
  ind.flex1 <- c(ind.flex1, ind.flex1.1)
  ind.flex3 <- c(ind.flex3, ind.flex3.1)
  ind.opt <- c(ind.opt, ind.opt.1)
  ind.jac <- c(ind.jac, ind.jac.1)
  ind.sim <- c(ind.sim, ind.sim.1)
  ind.ward <- c(ind.ward, ind.ward.1)
  ind.diana <- c(ind.diana, ind.diana.1)
  ind.kmeans <- c(ind.kmeans, ind.kmeans.1)
  ind.kulc <- c(ind.kulc, ind.kulc.1)
  ind.wardkulc <- c(ind.wardkulc, ind.wardkulc.1)
}
ind.table2 <- as.data.frame(cbind(klevel,ind.upgma,ind.flex,ind.flex1,ind.flex3,ind.opt,ind.jac,ind.sim,ind.ward,ind.diana,ind.kmeans,ind.kulc,ind.wardkulc))
ind.table2 <- ind.table2[-1,]
ind.table2 <- ind.table2 %>% mutate()
dni.table2 <- ind.table2[,-1] %>% t() %>% as.data.frame()
dni.table2 <- dni.table2 %>% mutate(s1to8 = apply(dni.table2[,1:8], MARGIN=1, FUN = 'mean'))
dni.table2 <- dni.table2 %>% mutate(s9to16 = apply(dni.table2[,9:16], MARGIN=1, FUN = 'mean'))





plotdata.total <- apply(plotdata2, MARGIN = 2, FUN = 'sum')
plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]



distbray <- vegdist(plotdata1, method='bray', binary=FALSE, na.rm=T)
distjac <- vegdist(plotdata1, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata1,method='Simpson'))
distkulc <- vegdist(plotdata1, method='kulczynski', binary=FALSE, na.rm=T)
tbrayagnes <- distbray %>% agnes(method = 'average')
tbrayflex05 <- distbray %>% flexbeta(beta= -0.05)
tbrayflex10 <- distbray %>% flexbeta(beta= -0.10)
tbrayflex15 <- distbray %>% flexbeta(beta= -0.15)
tbrayflex20 <- distbray %>% flexbeta(beta= -0.20)
tbrayflex25 <- distbray %>% flexbeta(beta= -0.25)
tbrayflex30 <- distbray %>% flexbeta(beta= -0.30)
tbrayflex35 <- distbray %>% flexbeta(beta= -0.35)
tbrayward <- distbray %>% agnes(method = 'ward')
k <- 2
klevel <- 0
sil.upgma <- 0
sil.flex05 <- 0
sil.flex10 <- 0
sil.flex15 <- 0
sil.flex20 <- 0
sil.flex25 <- 0
sil.flex30 <- 0
sil.flex35 <- 0
sil.ward <- 0

for (k in 2:20){
  sil.upgma.1 <- (tbrayagnes %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex05.1 <- (tbrayflex05 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex10.1 <- (tbrayflex10 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex15.1 <- (tbrayflex15 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex20.1 <- (tbrayflex20 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex25.1 <- (tbrayflex25 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex30.1 <- (tbrayflex30 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.flex35.1 <- (tbrayflex35 %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.ward.1 <- (tbrayward %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  
  klevel <- c(klevel, k)
  sil.upgma <- c(sil.upgma, sil.upgma.1)
  sil.flex05 <- c(sil.flex05, sil.flex05.1)
  sil.flex10 <- c(sil.flex10, sil.flex10.1)
  sil.flex15 <- c(sil.flex15, sil.flex15.1)
  sil.flex20 <- c(sil.flex20, sil.flex20.1)
  sil.flex25 <- c(sil.flex25, sil.flex25.1)
  sil.flex30 <- c(sil.flex30, sil.flex30.1)
  sil.flex35 <- c(sil.flex35, sil.flex35.1)
  sil.ward <- c(sil.ward, sil.ward.1)
}
sil.table <- as.data.frame(cbind(klevel,sil.upgma,sil.flex05,sil.flex10,sil.flex15,sil.flex20,sil.flex25,sil.flex30,sil.flex35,sil.ward))
sil.table <- sil.table[-1,]
sil.table <- sil.table %>% mutate()
lis.table <- sil.table[,-1] %>% t() %>% as.data.frame()
lis.table <- lis.table %>% mutate(s2to8 = apply(lis.table[,2:7], MARGIN=1, FUN = 'mean'))
lis.table <- lis.table %>% mutate(s9to16 = apply(lis.table[,8:16], MARGIN=1, FUN = 'mean'))

#----
plotdata.total <- apply(plotdata2, MARGIN = 2, FUN = 'sum')
plotdata.total <- as.data.frame(cbind(name=names(plotdata.total),total=plotdata.total))
removetaxon <- plotdata.total[plotdata.total$total %in% 0,]$name
plotdata1 <- plotdata[,!colnames(plotdata) %in% removetaxon]

distbray <- vegdist(plotdata1, method='bray', binary=FALSE, na.rm=T)# dbray <- as.data.frame(as.matrix(distbray))
distjac <- vegdist(plotdata1, method='jaccard', binary=FALSE, na.rm=T)
distsim <- as.dist(simil(plotdata1,method='Simpson'))
distkulc <- vegdist(plotdata1, method='kulczynski', binary=FALSE, na.rm=T)
tbrayagnes <- distbray %>% agnes(method = 'average')
tbrayflex05 <- distbray %>% flexbeta(beta= -0.05)
tbrayflex10 <- distbray %>% flexbeta(beta= -0.10)
tbrayflex15 <- distbray %>% flexbeta(beta= -0.15)
tbrayflex20 <- distbray %>% flexbeta(beta= -0.20)
tbrayflex25 <- distbray %>% flexbeta(beta= -0.25)
tbrayflex30 <- distbray %>% flexbeta(beta= -0.30)
tbrayflex35 <- distbray %>% flexbeta(beta= -0.35)
tbrayward <- distbray %>% agnes(method = 'ward')
k <- 2
klevel <- 0
ind.upgma <- 0
ind.flex05 <- 0
ind.flex10 <- 0
ind.flex15 <- 0
ind.flex20 <- 0
ind.flex25 <- 0
ind.flex30 <- 0
ind.flex35 <- 0
ind.ward <- 0

clu.upgma <- 0
clu.flex05 <- 0
clu.flex10 <- 0
clu.flex15 <- 0
clu.flex20 <- 0
clu.flex25 <- 0
clu.flex30 <- 0
clu.flex35 <- 0
clu.ward <- 0

timeA = Sys.time()
ind.upgma.1 <- optimclass(plotdata1, stride(8, as.hclust(tbrayagnes)))$sig.spc
Sys.time() - timeA
timeA = Sys.time()
# ind.upgma.1 <- tbrayagnes %>% cutree(k=8) 
# ind.upgma.1 <- indval(plotdata1, ind.upgma.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
# Sys.time() - timeA  
timeA = Sys.time()
for (k in 2:10){

  # ind.upgma.1 <- tbrayagnes %>% cutree(k=k) 
  # ind.flex05.1 <- tbrayflex05 %>% cutree(k=k)
  # ind.flex10.1 <- tbrayflex10 %>% cutree(k=k)
  # ind.flex15.1 <- tbrayflex15 %>% cutree(k=k)
  # ind.flex20.1 <- tbrayflex20 %>% cutree(k=k)
  # ind.flex25.1 <- tbrayflex25 %>% cutree(k=k)
  # ind.flex30.1 <- tbrayflex30 %>% cutree(k=k)
  # ind.flex35.1 <- tbrayflex35 %>% cutree(k=k)
  # ind.ward.1 <- tbrayward %>% cutree(k=k)
  # 
  # ind.upgma.1 <- indval(plotdata1, ind.upgma.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex05.1 <- indval(plotdata1, ind.flex05.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex10.1 <- indval(plotdata1, ind.flex10.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex15.1 <- indval(plotdata1, ind.flex15.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex20.1 <- indval(plotdata1, ind.flex20.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex25.1 <- indval(plotdata1, ind.flex25.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex30.1 <- indval(plotdata1, ind.flex30.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.flex35.1 <- indval(plotdata1, ind.flex35.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  # ind.ward.1 <- indval(plotdata1, ind.ward.1) %>% .$indval %>% apply(MARGIN=1, FUN = 'mean') %>% mean()
  
  ind.upgma.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayagnes)))
  ind.flex05.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex05)))
  ind.flex10.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex10)))
  ind.flex15.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex15)))
  ind.flex20.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex20)))
  ind.flex25.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex25)))
  ind.flex30.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex30)))
  ind.flex35.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayflex35)))
  ind.ward.0 <- optimclass(plotdata1, stride(k, as.hclust(tbrayward)))
  
  ind.upgma.1 <- ind.upgma.0$sig.spc
  ind.flex05.1 <- ind.flex05.0$sig.spc
  ind.flex10.1 <- ind.flex10.0$sig.spc
  ind.flex15.1 <- ind.flex15.0$sig.spc
  ind.flex20.1 <- ind.flex20.0$sig.spc
  ind.flex25.1 <- ind.flex25.0$sig.spc
  ind.flex30.1 <- ind.flex30.0$sig.spc
  ind.flex35.1 <- ind.flex35.0$sig.spc
  ind.ward.1 <- ind.ward.0$sig.spc
  
  clu.upgma.1 <- ind.upgma.0$sig.clust
  clu.flex05.1 <- ind.flex05.0$sig.clust
  clu.flex10.1 <- ind.flex10.0$sig.clust
  clu.flex15.1 <- ind.flex15.0$sig.clust
  clu.flex20.1 <- ind.flex20.0$sig.clust
  clu.flex25.1 <- ind.flex25.0$sig.clust
  clu.flex30.1 <- ind.flex30.0$sig.clust
  clu.flex35.1 <- ind.flex35.0$sig.clust
  clu.ward.1 <- ind.ward.0$sig.clust
  
  klevel <- c(klevel, k)
  ind.upgma <- c(ind.upgma, ind.upgma.1)
  ind.flex05 <- c(ind.flex05, ind.flex05.1)
  ind.flex10 <- c(ind.flex10, ind.flex10.1)
  ind.flex15 <- c(ind.flex15, ind.flex15.1)
  ind.flex20 <- c(ind.flex20, ind.flex20.1)
  ind.flex25 <- c(ind.flex25, ind.flex25.1)
  ind.flex30 <- c(ind.flex30, ind.flex30.1)
  ind.flex35 <- c(ind.flex35, ind.flex35.1)
  ind.ward <- c(ind.ward, ind.ward.1)
  
  clu.upgma <- c(clu.upgma, clu.upgma.1)
  clu.flex05 <- c(clu.flex05, clu.flex05.1)
  clu.flex10 <- c(clu.flex10, clu.flex10.1)
  clu.flex15 <- c(clu.flex15, clu.flex15.1)
  clu.flex20 <- c(clu.flex20, clu.flex20.1)
  clu.flex25 <- c(clu.flex25, clu.flex25.1)
  clu.flex30 <- c(clu.flex30, clu.flex30.1)
  clu.flex35 <- c(clu.flex35, clu.flex35.1)
  clu.ward <- c(clu.ward, clu.ward.1)
}
Sys.time() - timeA  
ind.table <- as.data.frame(cbind(klevel,ind.upgma,ind.flex05,ind.flex10,ind.flex15,ind.flex20,ind.flex25,ind.flex30,ind.flex35,ind.ward))
ind.table <- ind.table[-1,]
ind.table <- ind.table %>% mutate()
dni.table <- ind.table[,-1] %>% t() %>% as.data.frame()
dni.table <- dni.table %>% mutate(s2to8 = apply(dni.table[,1:7], MARGIN=1, FUN = 'mean'))
clu.table <- as.data.frame(cbind(klevel,clu.upgma,clu.flex05,clu.flex10,clu.flex15,clu.flex20,clu.flex25,clu.flex30,clu.flex35,clu.ward))
clu.table <- clu.table[-1,]
clu.table <- clu.table %>% mutate()
ulc.table <- clu.table[,-1] %>% t() %>% as.data.frame()
ulc.table <- ulc.table %>% mutate(s2to8 = apply(dni.table[,1:7], MARGIN=1, FUN = 'mean'))
write.csv(dni.table, 'output/optimclassindicators.csv', row.names = F)
write.csv(ulc.table, 'output/optimclasscluster.csv', row.names = F)
# dni.table2 <- dni.table
# write.csv(dni.table2, 'output/indval.csv', row.names = F)


ind <- indval(plotdata1, tbrayflex20 %>% cutree(k=4))
indicatorspecies <- as.data.frame(ind$relabu)