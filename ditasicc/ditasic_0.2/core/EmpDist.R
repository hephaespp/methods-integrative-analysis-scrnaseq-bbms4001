# Copyright (c) 2016,
# Martina Fischer, fischerm@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.



####   Pre-Step: Differential Testing Functions


# function < create.empdist > : Build empirical distribution based on poisson draws  ~ Taxa abundance distribution
# Input: 'abundance_objects' (output list of AbundanceEstimation.R) 
#        - absolute coefficient estimates  (as given by model, not normalized)
#        - absolute std.errors
#        - filtered taxa Ids  (mod.exclude)
#        - norm.fac: N2/N1  normalize for different number of input reads


create.empdist <- function(abundance_objects, norm.fac=1) {
  abs.coeff   <- abundance_objects[[1]]
  abs.error   <- abundance_objects[[2]]
  mod.exclude <- abundance_objects[[3]]
  lambdas <- list()
  filter <- ifelse(c(1:length(abs.coeff)) %in% mod.exclude, "FALSE", "TRUE")
  
  a.coeff <- abs.coeff[filter==T]/norm.fac          # estimates of abundant taxa
  a.error <- abs.error[filter==T]/sqrt(norm.fac)
  
  lower <- a.coeff - a.error
  upper <- a.coeff + a.error
  
  # 1.Step: Sample potential poisson lambdas from interval [estimate +/- std.error] for every taxa  
  set.seed(1442)
  N.lambda <- 100
  for (k in 1:length(a.coeff)) {
    if (lower[k] < 0) lower[k] <- 0 
    lambdas[[k]] <- runif(N.lambda, min = lower[k], max = upper[k]) 
  }
  
  # 2.Step:
  # for each abundant taxa (list entry): and for each of sampled lambda's (=mean read number of taxa): 
  # 500(N.poi) draws from lambda-defined poisson distribution -> specific lambda empirical distribution
  # Pooling all empiricial distributions of all lambdas -> overall empirical ditribution for each taxa
  set.seed(1442)
  N.poi <- 500
  empdist <- lapply(lambdas, function(i){as.vector(sapply(i, function(lambd) { rpois(N.poi, lambd) }))})
  
  
  # All taxa data: 'empdist.all' list object (0 = not ex. Taxa, 50000 poisson draws for existent taxa) 
  empdist.all <- sapply(names(abs.coeff), function(x) NULL)
  empdist.all[filter==F] <- 0
  empdist.all[filter==T] <- empdist
  
  return(empdist.all)
}
