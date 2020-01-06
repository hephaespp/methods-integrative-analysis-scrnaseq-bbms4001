countdc1_clred <- clr(countdc1)
#Problem: 'rmult' object created. How should it be converted back to matrix? 
#Action: Ignore the rmult created and directly work with the new object created 
hc <- hclust(dist(countdc1_clred), method = "complete") # data are clustered
#Problem: 'hierarchial clustering performed. How should the genes that belong to a certain cluster
#be recovered? 

#Temporary Solution:

params.groups <- newSplatParams(batchCells = 200, nGenes=1000) #parameters for control set
params.groups2 <- newSplatParams(batchCells = 200, nGenes=1000) #parameters for exp. cond. 1 set 
simdc1 <- splatSimulateGroups(params.groups, group.prob=c(0.25,0.25,0.5)) #simulation data for control set
simdc2 <- splatSimulateGroups(params.groups2, group.prob=c(0.15,0.15,0.7)) #simulation for the exp set 
countdc1 <- counts(simdc1)
pca <- prcomp(t(countdc1norm), rank = 2)

ggplot(as.data.frame(pca$x), aes(x=PC1, y=PC2)) + 
  geom_point() + 
  theme_bw() 

pcamat <- pca$x

km <- kmeans(pcamat, 3, 1000)
clust <- km$cluster
clust <- as.matrix(clust) 
select1 <- clust==1
select2 <- clust==2
select3 <- clust==3
clust1gene <- pcamat[select1, ]
clust2gene <- pcamat[select2, ]
clust3gene <- pcamat[select3, ]

