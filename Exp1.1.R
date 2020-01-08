# Load package
library(splatter)
library(scater)

####### Data Creation and Pre-processing  #######
params.groups <- newSplatParams(batchCells = 20000, nGenes=1000)
simdc <- splatSimulateGroups(params.groups, group.prob=c(0.25,0.25,0.25,0.25))
countdc <- counts(simdc)
normedsimdc <- normalize(simdc)
plotPCA(normedsimdc)

##Test2
install.packages("compositions")
library(compositions)
ilredcount <- ilr(t(countdc))
pca <- prcomp(ilredcount,rank=2)
ggplot(as.data.frame(pca$x),aes(x=PC1,y=PC2))+geom_point()+theme_bw()

normedcount <- sweep(countdc, 2, colSums(countdc), FUN="/")
pca1 <- prcomp(t(normedcount),rank=2)
ggplot(as.data.frame(pca1$x),aes(x=PC1,y=PC2))+geom_point()+theme_bw()

dfcount <- as.data.frame(countdc)
newdf <- sample_n(as.data.frame(t(dfcount)),2000)
newdf <- t(newdf)
normeddf <- sweep(newdf, 2, colSums(newdf), FUN="/")
pca2 <- prcomp(t(normeddf),rank=2)
ggplot(as.data.frame(pca2$x),aes(x=PC1,y=PC2))+geom_point()+theme_bw()

pcamat <- pca1$x
km <- kmeans(pcamat, 4, 1000)
clust <- km$cluster
clust <- as.matrix(clust) 
select1 <- clust==1
select2 <- clust==2
select3 <- clust==3
select4 <- clust==4

clust1gene <- pcamat[select1, ]
clust2gene <- pcamat[select2, ]
clust3gene <- pcamat[select3, ]
clust4gene <- pcamat[select4, ]

b <- list(nrow(clust1gene)/20000,nrow(clust2gene)/20000,nrow(clust3gene)/20000,nrow(clust4gene)/20000)
b <- as.numeric(b)
c <- b - c(0.25,0.25,0.25,0.25)
c <- abs(c)/0.25 *100
max(c)

