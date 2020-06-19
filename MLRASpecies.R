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
#ssecSpecies <- read.delim("data/USFSSubSpecies.txt")
#MLRASpecies <- read.delim("data/WardSpecies.txt")
#mlramatrix <- readRDS('data/mlramatrix.RDS')
#mlramatrix <- readRDS('data/usfsmatrix.RDS')
#fipsmatrix <- readRDS("data/fipsmatrix.RDS")
#ssecSpecies <- readRDS('data/USFSSubSpecies.RDS')
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
jactree <- agnes(fipsjacdist, method='average')
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
Species_km <- read.delim("data/Species_km.txt")
ssectmatrix <- readRDS('data/usfssubsecmatrix.RDS')
x1 <- which(Species_km$SumOfkm2 >= 1000000)
x2 <- which(Species_km$SumOfkm2 >= 200000 & Species_km$SumOfkm2 < 1000000)
x3 <- which(Species_km$SumOfkm2 < 200000)
jacdis1 <- as.data.frame(as.matrix(vegdist(ssectmatrix[,x1], method = 'jaccard')))
jacdis2 <- as.data.frame(as.matrix(vegdist(ssectmatrix[,x2], method = 'jaccard')))
jacdis3 <- as.data.frame(as.matrix(vegdist(ssectmatrix[,x3], method = 'jaccard')))
jactree1 <- agnes(jacdis1, method = 'average')
jactree2 <- agnes(jacdis2, method = 'average')
jactree3 <- agnes(jacdis3, method = 'average')
length(x1)
length(x2)
length(x3)
x <- as.data.frame(colnames(ssectmatrix[,x1]))

dcamodel <- decorana(ssectmatrix)
sites <- dcamodel$rproj
species <- dcamodel$cproj
dcadis1 <- as.data.frame(as.matrix(vegdist(sites[,1], method = 'euclidean')))
dcadis2 <- as.data.frame(as.matrix(vegdist(sites[,1:2], method = 'euclidean')))
dcadis3 <- as.data.frame(as.matrix(vegdist(sites[,1:3], method = 'euclidean')))
dcadis4 <- as.data.frame(as.matrix(vegdist(sites[,1:4], method = 'euclidean')))
dcadis <- (dcadis1*1 +dcadis2*1 +dcadis3*1 +dcadis4*1)/4
dcatree <- agnes(dcadis, method = 'average')


treeofinterest <- jactree2
cutsjac2 <- cutree(treeofinterest, k=2)
cutsjac4 <- cutree(treeofinterest, k=4)
cutsjac8 <- cutree(treeofinterest, k=8)
cutsjac16 <- cutree(treeofinterest, k=16)
cutsjac <- cutsjac2*10000+cutsjac4*1000+cutsjac8*100+cutsjac16
jactreelist <- cbind(as.data.frame(cutsjac2),as.data.frame(cutsjac4),as.data.frame(cutsjac8),as.data.frame(cutsjac16),as.data.frame(cutsjac), row.names(treeofinterest$data))
colnames(jactreelist) <- c('groups2', 'groups4','groups8','groups16','group','link')
write.dbf(jactreelist, 'output/subsectreelist.dbf')

sitesdf <- as.data.frame(sites)
sitesdf$link <- rownames(sitesdf)
write.dbf(sitesdf, 'output/dca.dbf')
w <- 500
h <- 10000
u <- 10
png(filename="output/subsecjac.png",width = w, height = h, units = "px", pointsize = u)
plot(as.phylo(as.hclust(jactree)), main='USFS subsection floristic simularity - jaccard metric',label.offset=0.05, direction='right', font=1, cex=0.85)
dev.off()

secdiversity <- cbind(rownames(ssectmatrix),(apply(ssectmatrix, MARGIN=1, FUN='sum')/100))
colnames(secdiversity) <- c('link', 'spp')
secdiversity <- aggregate(ssecSpecies[,c('abundance')], by=list(ssecSpecies$SUBSECTION), FUN='sum')
colnames(secdiversity) <- c('link', 'spp')
secdiversity$spp <- secdiversity$spp/100
secdiversity$link <- as.character(secdiversity$link)
secdiversity[substr(secdiversity$link,1,1) != 'M',]$link <- paste0('X',secdiversity[substr(secdiversity$link,1,1) != 'M',]$link)
write.dbf(secdiversity, 'output/secdiversity.dbf')

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




plotdata <- readRDS('data/fipsmatrix.RDS')

pcount <- apply(plotdata, MARGIN = 1, FUN = 'sum')

pcount2 <- pcount[pcount < 200 | grepl('F02',names(pcount))| grepl('F15',names(pcount))]
excluded <- names(pcount2)
#d <- vegdist(plotdata, method = 'bray')
#saveRDS(d, "data/fipsbraydist.RDS")
#d <- as.dist(simil(plotdata,method='Simpson'))
#saveRDS(d, "data/fipssimdist.RDS")

distbray <- readRDS("data/fipsbraydist.RDS")
distjac <- readRDS("data/fipsjacdist.RDS")
distsim <-  readRDS("data/fipssimdist.RDS")

plotdata <- plotdata[!rownames(plotdata) %in% excluded, !colnames(plotdata) %in% excluded]

distbray <- as.dist(as.matrix(distbray)[!rownames(as.matrix(distbray)) %in% excluded, !colnames(as.matrix(distbray)) %in% excluded])

distjac <- as.dist(as.matrix(distjac)[!rownames(as.matrix(distjac)) %in% excluded, !colnames(as.matrix(distjac)) %in% excluded])

distsim <- as.dist(as.matrix(distsim)[!rownames(as.matrix(distsim)) %in% excluded, !colnames(as.matrix(distsim)) %in% excluded])

tb <- distbray %>% agnes(method = 'average')
tj <- distjac %>% agnes(method = 'average')
tw <- distbray %>% agnes(method = 'ward') 
ts <- distsim %>% agnes(method = 'average') 
td <- distbray %>% diana 

silanalysis2 <- function(input){
  #distbray <- vegdist(input, method='bray', binary=FALSE, na.rm=T)
  #distjac <- vegdist(input, method='jaccard', binary=FALSE, na.rm=T)
  #distsim <- as.dist(simil(input,method='Simpson'))
  klist <- c(2:24,64,256)
  maxcluster <- min(24, nrow(input)-1)
  k <- 2
  klevel <- 0
  sil.bray <- 0
  sil.jac <- 0
  sil.sim <- 0
  sil.ward <- 0
  sil.diana <- 0
  
  for (k in 1:25){
    sil.bray1 <- (tb %>% cutree(k=klist[k]) %>% silhouette(distbray))[,3]%>% mean#[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    #sil.bray1 <- (aggregate(sil.bray1[,2], by=list(sil.bray1[,1]), FUN='mean'))[,2] %>% mean()
    sil.jac1 <- (tj %>% cutree(k=klist[k]) %>% silhouette(distbray))[,3]%>% mean#[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    #sil.jac1 <- (aggregate(sil.jac1[,2], by=list(sil.jac1[,1]), FUN='mean'))[,2] %>% mean()
    sil.sim1 <- (ts %>% cutree(k=klist[k]) %>% silhouette(distbray))[,3]%>% mean#[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    #sil.sim1 <- (aggregate(sil.sim1[,2], by=list(sil.sim1[,1]), FUN='mean'))[,2] %>% mean()
    sil.ward1 <- (tw %>% cutree(k=klist[k]) %>% silhouette(distbray))[,3]%>% mean#[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    #sil.ward1 <- (aggregate(sil.ward1[,2], by=list(sil.ward1[,1]), FUN='mean'))[,2] %>% mean()
    sil.diana1 <- (td %>% cutree(k=klist[k]) %>% silhouette(distbray))[,3]%>% mean#[,c(1,3)] %>% as.matrix()%>%as.data.frame() 
    #sil.diana1 <- (aggregate(sil.diana1[,2], by=list(sil.diana1[,1]), FUN='mean'))[,2] %>% mean()
    
    klevel <- c(klevel, klist[k])
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



k.agnes <- c(2,4,8,16,32,64,256)



agnes.cut1<- tj %>% cutree(k=k.agnes[1])
agnes.cut2<- tj %>% cutree(k=k.agnes[2])
agnes.cut3<- tj %>% cutree(k=k.agnes[3])
agnes.cut4<- tj %>% cutree(k=k.agnes[4])
agnes.cut5<- tj %>% cutree(k=k.agnes[5])
agnes.cut6<- tj %>% cutree(k=k.agnes[6])
agnes.cut7<- tj %>% cutree(k=k.agnes[7])

fips <- rownames(as.matrix(distbray))

fipsclust <- as.data.frame(cbind(fips,agnes.cut1,agnes.cut2,agnes.cut3,agnes.cut4,agnes.cut5,agnes.cut6,agnes.cut7))


write.dbf(fipsclust, 'output/fipsclust.dbf')

coph.agnes <- cor(cophenetic(as.hclust(ts)), cophenetic(as.hclust(tb)))
coph.jaccard <- cor(cophenetic(as.hclust(ts)), cophenetic(as.hclust(tj)))
coph.ward <- cor(cophenetic(as.hclust(ts)), cophenetic(as.hclust(tw)))
coph.diana <- cor(cophenetic(as.hclust(ts)), cophenetic(as.hclust(td)))
coph.sim <- cor(cophenetic(as.hclust(ts)), cophenetic(as.hclust(ts)))


coph.agnes
coph.jaccard
coph.ward
coph.diana
coph.sim

