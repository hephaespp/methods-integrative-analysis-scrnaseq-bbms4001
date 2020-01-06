####### Clustering the genes from coundc1 #######
clust1genevectors <- countdc1[,select1]
clust2genevectors <- countdc1[,select2]
clust3genevectors <- countdc1[,select3]
#create representative gene vector 
repre1 <- apply(clust1genevectors,1, FUN=median)
repre2 <- apply(clust2genevectors,1, FUN=median)
repre3 <- apply(clust3genevectors,1, FUN=median)
representatives <- cbind(repre1, repre2, repre3)
#pca to visualise distance
representativesnorm <- sweep(representatives, 2, colSums(representatives), FUN="/")
pca1 <- prcomp(t(representativesnorm),rank=2)
pcamat1 <- pca1$x
ggplot(as.data.frame(pca1$x),aes(x=PC1,y=PC2))+geom_point()+theme_bw()

#Row-variannce function 
RowVar <- function(x, ...) {
  rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
}
#Standard-deviations of each row in clust(i)genevectors 
sdclust1 <- sqrt(RowVar(clust1genevectors))
sdclust2 <- sqrt(RowVar(clust2genevectors))
sdclust3 <- sqrt(RowVar(clust3genevectors))
sdclust1 <- round(sdclust1)
sdclust2 <- round(sdclust2)
sdclust3 <- round(sdclust3)
#generate random gene vector characteristic of cluster 1 
generateclust1gene <- function() {
reprevector <- c(1:1000)
for (i in 1:1000) {
  delta <- sdclust1[i]
  e <- rpois(1,delta)
  reprevector[i] <- repre1[i] + e 
}
return(reprevector)
}

k <- list()
produceclust1gene <- function(k) {
Collection <- list()
  for (i in 1:k) {
    Collection[[i]] <- generateclust1gene()
  }
return(Collection)
}

bindgenestoprofiles <- function(G) {
a <- G[[1]]
  for (i in 2:length(G)) {
    a <- cbind(a,G[[i]])
  }
return(a)
}

X = c(8500,8300,8480,7960,8030)
Y = c(7710,7890,7920,8270.7860)
varx = var(X)
vary = var(Y)
r = sqrt((5*varx+5*vary)/8)*sqrt(2/5)
aa <- mean(X)-mean(Y)


