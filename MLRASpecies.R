# load required libraries
library(BiodiversityR)

library(plyr)
library(foreign)
library(vegan)
library(cluster)
library(ape)
######################################----

#MLRASpecies <- read.delim("data/MLRASpecies.txt")
mlramatrix <- readRDS('data/mlramatrix.RDS')
MLRASpecies <- readRDS('data/MLRASpecies.RDS')
#Matrix
#MLRASpecies$MLRA <-  as.character(paste0('M',MLRASpecies$MLRARSYM))
#mlramatrix <- makecommunitydataset(MLRASpecies, row = 'MLRA', column = 'Species', value = 'abundance')
#saveRDS(mlramatrix, "data/mlramatrix.RDS")
#saveRDS(MLRASpecies, "data/MLRASpecies.RDS")

mlrabraydist <- vegdist(mlramatrix,method='bray', binary=FALSE, na.rm=T)
mlrajacdist <- vegdist(mlramatrix,method='jaccard', binary=FALSE, na.rm=T)
mlrajacdist2 <- vegdist(mlramatrix,method='jaccard', binary=TRUE, na.rm=T)
mlrasorendist <- vegdist(mlramatrix,method='bray', binary=TRUE, na.rm=T)
disbray <- as.data.frame(as.matrix(mlrabraydist))
disjac <- as.data.frame(as.matrix(mlrajacdist))
disjac2 <- as.data.frame(as.matrix(mlrajacdist2))
dissoren <- as.data.frame(as.matrix(mlrasorendist))
jactree <- agnes(disjac, method='average')
jactree2 <- agnes(disjac2, method='average')
braytree <- agnes(disbray, method='average')
sorentree <- agnes(dissoren, method='average')
dianabraytree <- diana(disbray)
dianajactree <- diana(disbray)

treeofinterest <- jactree
cutsjac <- cutree(treeofinterest, k=15)
jactreelist <- cbind(as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('group','MLRA')
write.dbf(jactreelist, 'output/jacktreelist.dbf')


png(filename="output/jactree.png",width = 500, height = 2000, units = "px", pointsize = 12)

plot(as.phylo(as.hclust(jactree)), main='MLRA floristic simularity - jaccard',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/braytree.png",width = 500, height = 2000, units = "px", pointsize = 12)

plot(as.phylo(as.hclust(braytree)), main='MLRA floristic simularity - bray-curtis',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/dianajactree.png",width = 500, height = 2000, units = "px", pointsize = 12)

plot(as.phylo(as.hclust(dianajactree)), main='MLRA floristic simularity - diana jaccard',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/dianabraytree.png",width = 500, height = 2000, units = "px", pointsize = 12)

plot(as.phylo(as.hclust(dianabraytree)), main='MLRA floristic simularity - diana bray',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()


