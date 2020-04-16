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
######################################----

#subsections
plotdata <- readRDS('data/usfssubsecmatrix.RDS')
distbray <- readRDS('data/usfssubsec.distbray.RDS')
distjac <- readRDS('data/usfssubsec.distjac.RDS')
distsim <- readRDS('data/usfssubsec.distsim.RDS')


#validity using silhouette index
#----
#distbray <- vegdist(plotdata, method='bray', binary=FALSE, na.rm=T)
#distjac <- vegdist(plotdata, method='jaccard', binary=FALSE, na.rm=T)
#distsim <- as.dist(simil(plotdata,method='Simpson'))
#saveRDS(distbray,'data/usfssubsec.distbray.RDS')
#saveRDS(distjac,'data/usfssubsec.distjac.RDS')
#saveRDS(distsim,'data/usfssubsec.distsim.RDS')

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
for (k in 2:24){
  sil.bray1 <- (distbray %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.jac1 <- (distjac %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.sim1 <- (distsim %>% agnes(method = 'average') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.ward1 <- (distbray %>% agnes(method = 'ward') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.diana1 <- (distbray %>% diana %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.kmeans1 <- (kmeans(distbray, centers = k)$cluster %>% silhouette(distbray))[,3] %>% mean
  sil.single1 <- (distbray %>% agnes(method = 'single') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean
  sil.complete1 <- (distbray %>% agnes(method = 'complete') %>% cutree(k=k) %>% silhouette(distbray))[,3]%>% mean

  klevel <- c(klevel, k)
  sil.bray <- c(sil.bray, sil.bray1)
  sil.jac <- c(sil.jac, sil.jac1)
  sil.sim <- c(sil.sim, sil.sim1)
  sil.ward <- c(sil.ward, sil.ward1)
  sil.diana <- c(sil.diana, sil.diana1)
  sil.kmeans <- c(sil.kmeans, sil.kmeans1)
  sil.single <- c(sil.single, sil.single1)
  sil.complete <- c(sil.complete, sil.complete1)
}
sil.table <- as.data.frame(cbind(klevel,sil.bray,sil.jac,sil.sim,sil.ward,sil.diana,sil.kmeans,sil.single,sil.complete))
sil.table <- sil.table[-1,]
#saveRDS(sil.table,'output/usfssubsec.sil.table.RDS')


k.agnes <- c(2,3,4,7,10,16,24)
k.ward <- c(2,3,4,6,11,16,24)
k.diana <- c(2,3,4,7,10,16,24)
i<-1
agnes.cut <- 0
ward.cut <- 0
diana.cut <- 0

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


usfsclust <- as.data.frame(cbind(agnes.cut1,agnes.cut2,agnes.cut3,agnes.cut4,agnes.cut5,agnes.cut6,agnes.cut7,
                   ward.cut1,ward.cut2,ward.cut3,ward.cut4,ward.cut5,ward.cut6,ward.cut7,
                   diana.cut1,diana.cut2,diana.cut3,diana.cut4,diana.cut5,diana.cut6,diana.cut7))
usfsclust$subsect <- rownames(usfsclust)

write.dbf(usfsclust, 'output/usfsclust.dbf')
