# load required libraries
library(BiodiversityR)

library(plyr)
library(foreign)
library(vegan)
library(cluster)
library(ape)
library(proxy)
######################################----
MLRASpecies <- read.delim("data/USFSSubFamily.txt")
MLRASpecies <- read.delim("data/USFSSubHabitT.txt")
#MLRASpecies <- read.delim("data/EPA3Species.txt")
#MLRASpecies <- read.delim("data/EPASpecies.txt")
#MLRASpecies <- read.delim("data/MLRASpecies.txt")
#MLRASpecies <- read.delim("data/USFSSpecies.txt")
#MLRASpecies <- read.delim("data/USFSTrees.txt")
#MLRASpecies <- read.delim("data/USFSSubSpecies.txt")
#MLRASpecies <- read.delim("data/WardSpecies.txt")
#mlramatrix <- readRDS('data/mlramatrix.RDS')
#mlramatrix <- readRDS('data/usfsmatrix.RDS')
#fipsmatrix <- readRDS("data/fipsmatrix.RDS")
#MLRASpecies <- readRDS('data/USFSSubSpecies.RDS')
#MLRASpecies <- readRDS('data/WardSpecies.RDS')
#mlramatrix <- readRDS('data/usfsmatrix.RDS')
#mlramatrix <- readRDS('data/wardmatrix.RDS')
#fakeplots <- read.delim("data/fakeplots.txt")
#fipsplots <- read.delim("data/FIPSSpecies.txt")
#fipsplots <- readRDS("data/fipsplots.RDS")
mlramatrixfam <- mlramatrix
mlramatrixhab <- mlramatrix
jacdisthab <- jacdist <- as.data.frame(as.matrix(vegdist(mlramatrixhab,method='jaccard', binary=FALSE, na.rm=T)))
jacdistfam <- jacdist <- as.data.frame(as.matrix(vegdist(mlramatrixfam,method='jaccard', binary=FALSE, na.rm=T)))
jacdist2 <- jacdistfam+jacdisthab
#Matrix
mlramatrix <- makecommunitydataset(MLRASpecies, row = 'SUBSECTION', column = 'Habit', value = 'abun')
mlramatrix <- makecommunitydataset(MLRASpecies, row = 'SUBSECTION', column = 'Family', value = 'abundance')
#MLRASpecies$MLRA <-  as.character(paste0('M',MLRASpecies$MLRARSYM))
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'Level3', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'SUBSECTION', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'ward256', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'MLRA', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'USFS', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'SUBSECTION', column = 'Species', value = 'abundance')
#mlramatrix <- makecommunitydataset(fakeplots, row = 'plot', column = 'species', value = 'abundance')
#fipsmatrix <- makecommunitydataset(fipsplots, row = 'FFIPS', column = 'Species', value = 'abundance')
#saveRDS(fipsmatrix, "data/fipsmatrix.RDS")
#saveRDS(fipsplots, "data/fipsplots.RDS")
#saveRDS(mlramatrix, "data/epa3matrix.RDS")
#saveRDS(mlramatrix, "data/mlramatrix.RDS")
#saveRDS(mlramatrix, "data/treematrix.RDS")
#saveRDS(mlramatrix, "data/usfssubsecmatrix.RDS")
#saveRDS(mlramatrix, "data/wardmatrix.RDS")
#saveRDS(MLRASpecies, "data/MLRASpecies.RDS")
#saveRDS(MLRASpecies, "data/WardSpecies.RDS")
#saveRDS(MLRASpecies, "data/USFSSubSpecies.RDS")
#saveRDS(fipsjacdist, "data/fipsjacdist.RDS")
#fipsjacdist <- vegdist(fipsmatrix,method='jaccard', binary=FALSE, na.rm=T)
#subsecjacdist <- vegdist(mlramatrix,method='jaccard', binary=FALSE, na.rm=T)
#saveRDS(subsecjacdist, "data/subsecjacdist.RDS")

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
cutsjac4 <- cutree(treeofinterest, k=4)
cutsjac8 <- cutree(treeofinterest, k=8)
cutsjac16 <- cutree(treeofinterest, k=128)
cutsjac <- cutsjac2*10000+cutsjac4*1000+cutsjac8*100+cutsjac16
jactreelist <- cbind(as.data.frame(cutsjac2),as.data.frame(cutsjac4),as.data.frame(cutsjac8),as.data.frame(cutsjac16),as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('sgroups2', 'sgroups4','sgroups8','sgroups16','sgroup','MLRA')
#write.dbf(jactreelist, 'output/jacktreelist.dbf')
write.dbf(jactreelist, 'output/usfstreelist.dbf')

#FIPS
fipsjacdist <- readRDS("data/fipsjacdist.RDS")
fipsjacdistdf <- as.data.frame(as.matrix(fipsjacdist))
jactree <- agnes(fipsjacdistdf, method='average')
#saveRDS(jactree, 'data/fipsjactree.RDS')
treeofinterest <- jactree
cutsjac2 <- cutree(treeofinterest, k=2)
cutsjac4 <- cutree(treeofinterest, k=4)
cutsjac8 <- cutree(treeofinterest, k=8)
cutsjac16 <- cutree(treeofinterest, k=16)
cutsjac <- cutsjac2*10000+cutsjac4*1000+cutsjac8*100+cutsjac16
jactreelist <- cbind(as.data.frame(cutsjac2),as.data.frame(cutsjac4),as.data.frame(cutsjac8),as.data.frame(cutsjac16),as.data.frame(cutsjac), row.names(treeofinterest$data))
write.dbf(jactreelist, 'output/FIPStreelist.dbf')
w <- 500
h <- 20000
u <- 10
png(filename="output/fipsjac.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='fips floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()


#subsections

#mlramatrix <- readRDS('data/usfssubsecmatrix.RDS')
dcamodel <- decorana(mlramatrix)
sites <- dcamodel$rproj
species <- dcamodel$cproj
dcatree <- agnes(sites)
w <- 500
h <- 10000
u <- 10
png(filename="output/dcatree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(dcatree)), main='USFS subsection floristic simularity - DECORANA',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

#saveRDS(mlrasimpsim, "data/mlrasimpsim.RDS")
#mlrasimpsim <- readRDS("data/mlrasimpsim.RDS")
#subsecjacdist <- readRDS("data/subsecjacdist.RDS")
jacdist2 <- jacdistfam + jacdisthab + subsecjacdist
jacdist <- as.data.frame(as.matrix(vegdist(mlramatrix,method='jaccard', binary=FALSE, na.rm=T)))

jactree <- agnes(jacdist2, method='average')
#dianatree <- diana(MLRASpecies)
#simtree <- agnes(simsim, method='average')
#wardsimtree <- agnes(simsim, method='ward')
#kmeangroups <- kmeans(fipsjacdistdf, centers = 256)
#kclust <- as.data.frame( kmeangroups$cluster)
#colnames(kclust) <- 'clust'
#kclust$link <- rownames(kclust)
#write.dbf(kclust, 'output/kclust.dbf')
#saveRDS(jactree, 'data/fipsjactree.RDS')
treeofinterest <- jactree
cutsjac2 <- cutree(treeofinterest, k=32)
cutsjac4 <- cutree(treeofinterest, k=64)
cutsjac8 <- cutree(treeofinterest, k=128)
cutsjac16 <- cutree(treeofinterest, k=256)
cutsjac <- cutsjac2*10000+cutsjac4*1000+cutsjac8*100+cutsjac16
jactreelist <- cbind(as.data.frame(cutsjac2),as.data.frame(cutsjac4),as.data.frame(cutsjac8),as.data.frame(cutsjac16),as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('groups2', 'groups4','groups8','groups16','group','link')
write.dbf(jactreelist, 'output/subsectreelist.dbf')
w <- 500
h <- 10000
u <- 10
png(filename="output/subsecjac.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='USFS subsection floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()


#wardclusters
mlramatrix <- readRDS('data/wardmatrix.RDS')
mlrajacdist <- vegdist(mlramatrix,method='jaccard', binary=FALSE, na.rm=T)
mlrasimpsim <- simil(mlramatrix,method='Simpson')
simsim <- as.data.frame(as.matrix(mlrasimpsim))
wardjacdistdf <- as.data.frame(as.matrix(mlrajacdist))
jactree <- agnes(wardjacdistdf, method='average')
simtree <- agnes(simsim, method='average')
wardtree <- agnes(wardjacdistdf, method='ward')
treeofinterest <- jactree
cutsjac2 <- cutree(treeofinterest, k=2)
cutsjac4 <- cutree(treeofinterest, k=4)
cutsjac8 <- cutree(treeofinterest, k=8)
cutsjac16 <- cutree(treeofinterest, k=16)
cutsjac <- cutsjac2*10000+cutsjac4*1000+cutsjac8*100+cutsjac16
jactreelist <- cbind(as.data.frame(cutsjac2),as.data.frame(cutsjac4),as.data.frame(cutsjac8),as.data.frame(cutsjac16),as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('wgroups2', 'wgroups4','wgroups8','wgroups16','wgroup','wlink')
write.dbf(jactreelist, 'output/wardtreelist.dbf')
w <- 500
h <- 10000
u <- 10
png(filename="output/subsecjac.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='USFS subsection floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

#epa

#mlramatrix <- readRDS('data/epamatrix.RDS')
#epajacdist <- as.data.frame(as.matrix(vegdist(mlramatrix,method='jaccard', binary=FALSE, na.rm=T)))
#saveRDS(epajacdist,'data/epajacdist.RDS')
epajacdist <- readRDS("data/epajacdist.RDS")
#saveRDS(mlramatrix, 'data/subsecsimpson.RDS')
jactree <- agnes(epajacdist, method='average')
dianajactree <- diana(epajacdist)
#dianatree <- diana(MLRASpecies)
#simtree <- agnes(simsim, method='average')
treeofinterest <- jactree
cutsjac2 <- cutree(treeofinterest, k=32)
cutsjac4 <- cutree(treeofinterest, k=64)
cutsjac8 <- cutree(treeofinterest, k=128)
cutsjac16 <- cutree(treeofinterest, k=256)
cutsjac <- cutsjac2*10000+cutsjac4*1000+cutsjac8*100+cutsjac16
jactreelist <- cbind(as.data.frame(cutsjac2),as.data.frame(cutsjac4),as.data.frame(cutsjac8),as.data.frame(cutsjac16),as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('sgroups2', 'sgroups4','sgroups8','sgroups16','sgroup','link')
write.dbf(jactreelist, 'output/epa3treelist.dbf')
w <- 500
h <- 10000
u <- 10
png(filename="output/epajactree.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='EPA floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()


