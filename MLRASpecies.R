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

mlrabraydist <- vegdist(mlramatrix,method='bray', na.rm=T)
mlrajacdist <- vegdist(mlramatrix,method='jaccard', na.rm=T)
disbray <- as.data.frame(as.matrix(mlrabraydist))
disjac <- as.data.frame(as.matrix(mlrajacdist))
jactree <- agnes(disjac, method='average')
braytree <- agnes(disbray, method='average')
png(filename="output/jactree.png",width = 500, height = 2000, units = "px", pointsize = 12)

plot(as.phylo(as.hclust(jactree)), main='MLRA floristic simularity - jaccard',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

png(filename="output/braytree.png",width = 500, height = 2000, units = "px", pointsize = 12)

plot(as.phylo(as.hclust(braytree)), main='MLRA floristic simularity - bray-curtis',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

