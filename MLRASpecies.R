# load required libraries
library(BiodiversityR)

library(plyr)
library(foreign)
library(vegan)
library(cluster)
library(ape)
library(proxy)
######################################----

#MLRASpecies <- read.delim("data/MLRASpecies.txt")
#MLRASpecies <- read.delim("data/USFSSpecies.txt")
mlramatrix <- readRDS('data/mlramatrix.RDS')
#mlramatrix <- readRDS('data/usfsmatrix.RDS')
#MLRASpecies <- readRDS('data/MLRASpecies.RDS')
#MLRASpecies <- readRDS('data/usfsmatrix.RDS')
fakeplots <- read.delim("data/fakeplots.txt")
#fipsplots <- read.delim("data/FIPSSpecies.txt")

#Matrix
#MLRASpecies$MLRA <-  as.character(paste0('M',MLRASpecies$MLRARSYM))
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'MLRA', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'USFS', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(fakeplots, row = 'plot', column = 'species', value = 'abundance')
#fipsmatrix <- makecommunitydataset(fipsplots, row = 'FFIPS', column = 'Species', value = 'abundance')

#saveRDS(mlramatrix, "data/mlramatrix.RDS")
#saveRDS(mlramatrix, "data/usfsmatrix.RDS")
#saveRDS(MLRASpecies, "data/MLRASpecies.RDS")

mlrabraydist <- vegdist(mlramatrix,method='bray', binary=FALSE, na.rm=T)
mlrajacdist <- vegdist(mlramatrix,method='jaccard', binary=FALSE, na.rm=T)
mlrajacdist2 <- vegdist(mlramatrix,method='jaccard', binary=TRUE, na.rm=T)
mlrasorendist <- vegdist(mlramatrix,method='bray', binary=TRUE, na.rm=T)
mlrasimpsim <- simil(mlramatrix,method='Simpson')
disbray <- as.data.frame(as.matrix(mlrabraydist))
disjac <- as.data.frame(as.matrix(mlrajacdist))
disjacsqrt <- as.data.frame(as.matrix(vegdist(mlramatrix^(1/3),method='jaccard', binary=FALSE, na.rm=T)))
disbraysqrt <- as.data.frame(as.matrix(vegdist(mlramatrix^(1/3),method='bray', binary=FALSE, na.rm=T)))
disjac2 <- as.data.frame(as.matrix(mlrajacdist2))
dissoren <- as.data.frame(as.matrix(mlrasorendist))
simsim <- as.data.frame(as.matrix(mlrasimpsim))
chaodist <- as.data.frame(as.matrix(vegdist(mlramatrix,method='chao', binary=FALSE, na.rm=T)))
#tree build
jactree <- agnes(disjac, method='average')
jacsqrttree <- agnes(disjacsqrt, method='average')
jactree2 <- agnes(disjac2, method='average')
braytree <- agnes(disbray, method='average')
braysqrttree <- agnes(disbraysqrt, method='average')
sorentree <- agnes(dissoren, method='average')
dianabraytree <- diana(disbray)
dianajactree <- diana(disjac)
simtree <- agnes(simsim, method='average')
chaotree <- agnes(chaodist, method='average')
w <- 500
h <- 2000
u <- 12

png(filename="output/jactree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='MLRA floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()
png(filename="output/jacsroottree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jacsqrttree)), main='MLRA floristic simularity - jaccard cuberoot',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()
png(filename="output/jactree2.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree2)), main='MLRA floristic simularity - jaccard binary',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/braytree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(braytree)), main='MLRA floristic simularity - bray-curtis',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()
png(filename="output/brayroottree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(braysqrttree)), main='MLRA floristic simularity - bray-curtis-cuberoot',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()
png(filename="output/sorentree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(sorentree)), main='MLRA floristic simularity - sorensen',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/dianajactree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(dianajactree)), main='MLRA floristic simularity - diana jaccard',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/dianabraytree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(dianabraytree)), main='MLRA floristic simularity - diana bray',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/simsim.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(simtree)), main='MLRA floristic simularity - simpson',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()
png(filename="output/chaotree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(chaotree)), main='MLRA floristic simularity - chao',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()


treeofinterest <- jactree
cutsjac2 <- cutree(treeofinterest, k=2)
cutsjac5 <- cutree(treeofinterest, k=6)
cutsjac10 <- cutree(treeofinterest, k=15)
cutsjac <- cutsjac2*1000+cutsjac5*100+cutsjac10
jactreelist <- cbind(as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('group','MLRA')
write.dbf(jactreelist, 'output/jacktreelist.dbf')




