# Copyright (c) 2016,
# Martina Fischer, fischerm@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.


## Differential Abundance Testing 


# Function < DiffExpr > :  Difference distribution approach for empirical p-value calculation

# Input:  - empdist objects (output of EmpDist.R) for the two samples considered
#         - normalized abundance estimates of sample 1+2 
#         - threshold.minreads: minimum number of reads to assign significant taxa existence
#         - seed: default for internal sampling process


DiffExpr  <- function(empdist1, empdist2, norm.coeff_s1, norm.coeff_s2, threshold.minreads =0, seed=1448) {
  
  glm.coeff.1 <- norm.coeff_s1
  glm.coeff.2 <- norm.coeff_s2 
  N.draws <- 50000                                              # N.draws = N.lambda * N.poi (defined in <create.empdist>)
  
  result <- data.frame("norm.count.estimate.1"=NA, "norm.count.estimate.2"=NA, "diff.pval"=NA, "pval.adj"=NA, "FC.log2"=NA)
  
  for (k in 1:length(glm.coeff.1))   {
    
    if (length(empdist1[[k]])==1) 
    { result[k, "norm.count.estimate.1"] <- NA;
    if (length(empdist2[[k]])==1) {result[k, "norm.count.estimate.2"] <- NA}   # case: o o
    else {                                                      # case: o x
      result[k, "norm.count.estimate.2"] <- glm.coeff.2[k]                     # normalized glm estimates
      diff.emp <- empdist2[[k]]- threshold.minreads             # sig. for existence: > min.threshold input reads
      pval2 <- ecdf(diff.emp)(0)
      result[k, "diff.pval"] <- pval2
    }
    }
    else {
      if (length(empdist2[[k]])==1)                             # case: x o
      {result[k, "norm.count.estimate.2"] <- NA;
      result[k, "norm.count.estimate.1"] <- glm.coeff.1[k]                     # normalized glm estimates 
      diff.emp <- empdist1[[k]]- threshold.minreads             # sig. for existence: > min.threshold input reads
      pval1 <- ecdf(diff.emp)(0)
      result[k, "diff.pval"] <- pval1
      }  
      else {                                                    # case: x x 
        result[k, "norm.count.estimate.1"] <- glm.coeff.1[k]                   # normalized glm estimates
        result[k, "norm.count.estimate.2"] <- glm.coeff.2[k]
        result[k, "FC.log2"] <- log2(glm.coeff.2[k]/glm.coeff.1[k])
        
        set.seed(seed)
        emp1 <- sample(empdist1[[k]], N.draws)
        emp2 <- sample(empdist2[[k]], N.draws)
        diff.emp <- emp1 - emp2                                 # empirical difference distribution
        
        pval <- ecdf(diff.emp)(0)
        if (pval > 0.5) pval <- 1- pval
        result[k, "diff.pval"] <- pval
      }
    }
  } # end for-loop
    
  result[, "pval.adj"] <- p.adjust(result[, "diff.pval"], "BH")
  return(result)
}



# print results

print_output <- function(result, taxa, filename){
  
  result <- cbind("taxa.name"=as.vector(taxa), result)
  result[,c(2,3,6)] <- round(result[,c(2,3,6)],6)
  write.table(result, file=filename, quote=F, sep="\t", row.names=F)
}
