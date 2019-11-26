# load required libraries
library(BiodiversityR)

library(plyr)
library(foreign)
library(vegan)
library(cluster)
library(ape)
library(proxy)
######################################----
ssecspp <- readRDS('data/USFSSubSpecies.RDS')
ssectmatrix<- readRDS('data/usfssubsecmatrix.RDS')

ssecspp$SUBSECTION <- as.character(ssecspp$SUBSECTION)
ssecspp[substr(ssecspp$SUBSECTION,1,1) != 'M',]$SUBSECTION <- paste0('X', ssecspp[substr(ssecspp$SUBSECTION,1,1) != 'M',]$SUBSECTION)
resultsdf <- as.data.frame(cbind(
  methodsx=c('jactree','bratree','krutree', 'jwardree','bwardree','jdiatree','bdiatree'),
  kcut002=c(0,0,0,0,0,0,0),
  kcut004=c(0,0,0,0,0,0,0),
  kcut008=c(0,0,0,0,0,0,0),
  kcut016=c(0,0,0,0,0,0,0),
  kcut032=c(0,0,0,0,0,0,0),
  kcut064=c(0,0,0,0,0,0,0)))

for (i in 2:7){
resultsdf[,i] <- as.numeric(resultsdf[,i])
resultsdf[,i] <- NA}

  

jacdist <- as.data.frame(as.matrix(vegdist(ssectmatrix,method='jaccard', binary=FALSE, na.rm=T)))
bradist <- as.data.frame(as.matrix(vegdist(ssectmatrix,method='bray', binary=FALSE, na.rm=T)))
krudist <- as.data.frame(as.matrix(vegdist(ssectmatrix,method='kulczynski', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist, method = 'average')
bratree <- agnes(bradist, method = 'average')
krutree <- agnes(krudist, method = 'average')
jwardree <- agnes(jacdist, method = 'ward')
bwardree <- agnes(bradist, method = 'ward')
jdiatree <- diana(jacdist)
bdiatree <- diana(bradist)

treeofinterest <- list(jactree,bratree,krutree,jwardree,bwardree,jdiatree,bdiatree)
precutlist <- c(2,4,8,16,32,64)
for (i in 1:7){
for (j in 1:5){
cuts000 <- cutree(treeofinterest[[i]], k=precutlist[j])

cuttreelist <- cbind(as.data.frame(cuts000), row.names(treeofinterest[[i]]$data))
colnames(cuttreelist) <- c('cuts000','ssec')
ssecsppclust <- merge(ssecspp, cuttreelist, by.x = 'SUBSECTION', by.y = 'ssec')

ssecsppclust.mean <- aggregate(ssecsppclust[,"abundance"], by=list(ssecsppclust$Species, ssecsppclust$cuts000), FUN='sum')
colnames(ssecsppclust.mean) <- c('taxon', 'cluster', 'sum')
ssecsppclust.count <- aggregate(unique(ssecsppclust[c('cuts000', 'SUBSECTION')])$SUBSECTION, 
                                 by=list(unique(ssecsppclust[c('cuts000', 'SUBSECTION')])$cuts000), FUN='length')
colnames(ssecsppclust.count) <- c('cluster', 'count')
ssecsppclust.mean <- merge(ssecsppclust.mean, ssecsppclust.count, by='cluster')
ssecsppclust.mean$mean <- ssecsppclust.mean$sum/ssecsppclust.mean$count

ssecsppclust.mean.aff <- aggregate(ssecsppclust.mean[,"mean"], by=list(ssecsppclust.mean$taxon), FUN='sum')
colnames(ssecsppclust.mean.aff) <- c('taxon', 'affsum')
ssecsppclust.mean <- merge(ssecsppclust.mean, ssecsppclust.mean.aff, by = 'taxon')
ssecsppclust.mean$aff <- ssecsppclust.mean$mean/ssecsppclust.mean$affsum*100
ssecsppclust.mean.aff.max <- aggregate(ssecsppclust.mean[,"aff"], by=list(ssecsppclust.mean$taxon), FUN='max')
colnames(ssecsppclust.mean.aff.max) <- c('taxon', 'max')
meanaff <- mean(ssecsppclust.mean.aff.max$max)
resultsdf[i,j+1] <- meanaff
}}
#----
library(clusterCrit)
for (i in 1:7){
  for (j in 1:6){
    cuts000 <- cutree(treeofinterest[[i]], k=precutlist[j])
    #resultsdf[i,j+1] <- cor(as.dist(jacdist), cophenetic(treeofinterest[[i]]))
    
    resultsdf[i,j+1] <- intCriteria(as.matrix(bradist), cuts000, "Dunn")
  }}

cop_jaba <- cor(as.dist(jacdist), cophenetic(bratree))
cop_baja <- cor(as.dist(bradist), cophenetic(jactree))
cop_jac <- cor(as.dist(jacdist), cophenetic(jactree))
cop_bra <- cor(as.dist(bradist), cophenetic(bratree))
cop_kruc <- cor(as.dist(krudist), cophenetic(krutree))
cop_brward <- cor(as.dist(bradist), cophenetic(bwardree))
cop_jacward <- cor(as.dist(jacdist), cophenetic(jwardree))
cop_brdian <- cor(as.dist(bradist), cophenetic(bdiatree))
cop_jacdian <- cor(as.dist(jacdist), cophenetic(jdiatree))
#----
comparedf <- as.data.frame(cbind(
  methodsx=c('jactree','bratree','krutree', 'jwardree','bwardree','jdiatree','bdiatree'),
  kcut002=c(0,0,0,0,0,0,0),
  kcut008=c(0,0,0,0,0,0,0),
  kcut032=c(0,0,0,0,0,0,0)))

for (i in 2:7){
  resultsdf[,i] <- as.numeric(resultsdf[,i])
  resultsdf[,i] <- NA}

precutlist <- c(2,8,32)

jacdist.mean <- apply(jacdist, MARGIN = 1, FUN ='mean')
jacdist.mean <- mean(jacdist.mean)
for (i in 1:7){
  for (j in 1:3){
    cuts000 <- cutree(treeofinterest[[i]], k=precutlist[j])
    
    cuttreelist <- cbind(as.data.frame(cuts000), row.names(treeofinterest[[i]]$data))
    colnames(cuttreelist) <- c('cuts000','ssec')
    ssecsppclust <- merge(ssecspp, cuttreelist, by.x = 'SUBSECTION', by.y = 'ssec')
    
    ssecsppclust.mean <- aggregate(ssecsppclust[,"abundance"], by=list(ssecsppclust$Species, ssecsppclust$cuts000), FUN='sum')
    colnames(ssecsppclust.mean) <- c('taxon', 'cluster', 'sum')
    ssecsppclust.count <- aggregate(unique(ssecsppclust[c('cuts000', 'SUBSECTION')])$SUBSECTION, 
                                    by=list(unique(ssecsppclust[c('cuts000', 'SUBSECTION')])$cuts000), FUN='length')
    colnames(ssecsppclust.count) <- c('cluster', 'count')
    ssecsppclust.mean <- merge(ssecsppclust.mean, ssecsppclust.count, by='cluster')
    ssecsppclust.mean$mean <- ssecsppclust.mean$sum/ssecsppclust.mean$count
    listclust <- unique(ssecsppclust.mean$cluster)
    nclust <- length(listclust)
    prejac <- 0
    for (k in 1:nclust){
      for (m in 1:nclust){
        A <- ssecsppclust.mean[ssecsppclust.mean$cluster == listclust[k],c('taxon', 'mean')]
        B <- ssecsppclust.mean[ssecsppclust.mean$cluster == listclust[m],c('taxon', 'mean')]
        C <- merge(A,B, by='taxon', all.x = TRUE,all.y = TRUE)
        C$mean.x <- ifelse(is.na(C$mean.x), 0, C$mean.x)
        C$mean.y <- ifelse(is.na(C$mean.y), 0, C$mean.y)
        C$abdif <- abs(C$mean.x-C$mean.y)
        C <- sum(C$abdif)
        A <- sum(A$mean)
        B <- sum(B$mean)
        bray <- C/(A+B)
        jaccard <- 2*bray/(1+bray)
        prejac <- jaccard + prejac
      }}
    comparedf[i,j+1] <- prejac/nclust^2
  }}

#saveRDS(resultsdf,'output/resultdf.RDF')