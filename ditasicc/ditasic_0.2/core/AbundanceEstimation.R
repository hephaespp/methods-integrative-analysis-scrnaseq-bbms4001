# Copyright (c) 2016,
# Martina Fischer, fischerm@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.


####   Abundance Estimation


# Input required: - mat: path to 'similarity_matrix.npy'                     (Output ditasic_matrix.py)
#                 - mapped.read.counts: path to 'sample_count.npy'           (Output ditasic_mapping.py )
#                 - N (total number of input reads): sample_total.npy        (Output ditasic_mapping.py ) 
# optional arguments: filtering (logical), p-value-threshold, min.threshold 


AbundanceEstim <- function( mat, mapped.read.counts, N, filtering=F, pval.threshold=0.05, min.threshold=0 ) {
  # READ - IN
  # Similarity Matrix
  sim.mat <- mat
  sim.matn <- mat.normtranspose(sim.mat)        # transponation + normalization
  
  if(FALSE %in% sapply(1:dim(sim.matn)[1],function(i){sim.matn[i,i]==1}) ){
    stop("Error in matrix construction, entries in matrix diagonal are expected to be 1 after normalization")
  }    
  
  # observed mapping counts
  read.count <- mapped.read.counts
  N <- N
  

  # MODELLING STEP    
  # numerical stabilization
  zero.counts <- which(read.count==0)
  if (length(zero.counts) >0) {
    read.count[zero.counts] <- 1e-10
  }
  if (length(zero.counts) > length(read.count)/3) {
    warning(" More than a third of the observed mapping counts are zero. This likely reflects non-existence of many 'species' selected for this analysis and can cause an estimation bias. It is recommended to specify a set of expected species in the data set with more focus on its underlying strains of interest.")
  }
  
  # FULL Model Poisson + Identity Link
  glm.model <- glm(read.count ~ 0+ sim.matn, family="poisson"(link = "identity"), maxit=100, trace=T)
  
  
  # Overdispersion test
  if(dispersiontest(glm.model)$p.value < 0.05)  {  
  glm.model <- glm(read.count ~ 0+ sim.matn, family="quasipoisson"(link = "identity"), maxit=100, trace=T)
  }
  
  
  # estimates:
  glm.coeffs <- glm.model$coefficients
  error.abs <- summary(glm.model)$coefficients[,2]
  raw.pval <- summary(glm.model)$coefficients[,4]
  
  # Filtering by p-value and min.threshold
  if (filtering == T){
    mod.exclude1 <- which(summary(glm.model)$coefficients[,4] > pval.threshold) 
    mod.exclude2 <- which(abs(glm.coeffs) < min.threshold)  
    mod.exclude<- sort(union(mod.exclude1, mod.exclude2))
    glm.coeffs[mod.exclude] <- 0
  }
  else mod.exclude <- NULL
  

  # re-optimization (case of negative estimates)
  if(sum(glm.coeffs <0) >0)
  { glm.coeffs <- optim.fit(read.count, sim.matn, start=glm.coeffs)$par
  }
  
  return(list(glm.coeffs, error.abs, mod.exclude, raw.pval ))
} 
#function end




# Similiarity matrix normalization + transponation
mat.normtranspose <- function(mat){
  matn  <- matrix(NA, dim(mat)[1],dim(mat)[1])
  for (i in 1:dim(mat)[1]){
    for (j in 1:dim(mat)[2]){
      matn[i,j] <- mat[j,i]/mat[i,i]
    }
  }
  return(matn)
}


# Optimization function
optim.fit <- function(y, z, start=rep(mean(y), ncol(z)),...){    
  deviance <- function(beta, y, z){
    mu <- z %*% beta; 2*sum(mu - y - y*log(mu/y))
  }
  grad<- function(beta, y, z){
    mu <- z %*% beta; 2* t(1 - y/mu) %*% z   
  }
  nlminb(start, deviance, grad, lower =0, y=y, z=z,  ...) 
}


# print results

Abundance_output <- function(abund_result, taxa, filename){
  
  filter.vec <- rep("no",length(as.vector(taxa)))
  filter.vec[abund_result[[3]]] <- "yes"
  result <- data.frame("taxa.name"=as.vector(taxa), "count.estimate"=as.numeric(as.vector(abund_result[[1]])), "error.estimate"=as.numeric(as.vector(abund_result[[2]])), "filtered"=filter.vec, "raw.pval"=as.numeric(as.vector(abund_result[[4]])))
  result[,c(2,3)] <- round(result[,c(2,3)],2)
  result[,5] <- round(result[,5],6)
  write.table(result, file=filename, quote=F, sep="\t", row.names=F)
}
