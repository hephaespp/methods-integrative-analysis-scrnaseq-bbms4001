# Copyright (c) 2016,
# Martina Fischer, fischerm@rki.de, Robert Koch Institute, Germany, 
# All rights reserved. For details, please note the license.txt.



# Data load and file checks

# Read-In npy.files (outputs from 'ditasic_mapping' and 'distasic_matrix')


data.load <- function(refs, mat, counts_s1, N_s1, counts_s2, N_s2){
  
  # read-in taxa names
  taxa <- read.table(refs, col.names="V", row.names=NULL, colClasses = "character")$V
  taxa <- sapply(taxa, function(k) {a <- strsplit(k, "/")[[1]]; as.character(a[length(a)]) }, USE.NAMES=F)
  
  #similarity matrix
  sim.mat  <- npyLoad(mat)
  
  # observed mapping counts / input number - sample 1
  read.count.s1 <- npyLoad(counts_s1)                
  N.s1 <- npyLoad(N_s1)               
  
  # observed mapping counts / input number - sample 2
  if(is.na(counts_s2)!=TRUE && is.na(N_s2)!=TRUE ) 
  {
  read.count.s2 <- npyLoad(counts_s2)                
  N.s2 <- npyLoad(N_s2)       
  }
  else {
  read.count.s2 <- NA
  N.s2 <- NA
  }
  
  
  # Data check 
  # stop executions if requirements are not fulfilled
  # value check:
  if (all(sim.mat>=0) != TRUE) {
    stop("Error in similarity matrix, entries need to positive integer values (including 0)")
  }
  if (sign(N.s1) != 1 ) {
    stop("Error in total number of input reads (sample 1), value need to be a positive integer")
  }
  if (all(read.count.s1>=0) != TRUE ) {
    stop("Error in mapped count vector (sample 1), entries need to positive integer values")
  }
  
  if(is.na(counts_s2)!=TRUE && is.na(N_s2)!=TRUE ) 
  {
	if (sign(N.s2) != 1) {
    stop("Error in total number of input reads (sample 2), value need to be a positive integer") }
	if (all(read.count.s2>=0) != TRUE) {
    stop("Error in mapped count vector (sample 2), entries need to positive integer values")    }
	}
  
  
  # dimension check:
  if (length(read.count.s1) != dim(sim.mat)[1] ) {
    stop("'Number of observed mapping counts (sample 1) must equal the number of taxa used for construction of the similarity matrix ' ")
  }
  if(length(taxa) != dim(sim.mat)[1]) {
    stop(" Number of taxa in taxa-file does not equal number of taxa in matrix calculation ")
  }
  
  if(is.na(counts_s2)!=TRUE && is.na(N_s2)!=TRUE ) 
  {
	if (length(read.count.s2) != dim(sim.mat)[1]) {
    stop("'Number of observed mapping counts (sample 2) must equal the number of taxa used for construction of the similarity matrix ' ") }
	}
	
	
  # mapping check:
  if (sum(read.count.s1) < (N.s1*0.75)) {
    stop("'More than 25% of the reads could not be assigned, check mapping accuracy and suitability of selected reference genomes ': ",sum(read.count.s1), " < ",N.s1, " (total number of reads in Sample 1)")
  }
  
  if(is.na(counts_s2)!=TRUE && is.na(N_s2)!=TRUE ) 
  {
  if (sum(read.count.s2) < (N.s2*0.75)) {
    stop("'More than 25% of the reads could not be assigned, check mapping accuracy and suitability of selected reference genomes': ",sum(read.count.s2), " < ",N.s2, " (total number of reads in Sample 2)")   }
  }
  
  
  
  # resolution check:
  if ((sum(read.count.s1) < N.s1) && (sum(read.count.s1) > (N.s1*0.75))) {
    stop("'Sum of observed mapping counts (sample 1) does not exceed the total number of reads. Purpose of the method is the resolution of shared read counts. In case of purely unique mappings, taxa abundance estimates can be directly inferred from mapping counts. ': ",sum(read.count.s1), " < ",N.s1, " (total number of reads in Sample 1)")
  }
  
  if(is.na(counts_s2)!=TRUE && is.na(N_s2)!=TRUE ) 
  {
  if ((sum(read.count.s2) < N.s2) && (sum(read.count.s2) < (N.s2*0.75))) { 
    stop("'Sum of observed mapping counts (sample 2) does not exceed the total number of reads. Purpose of the method is the resolution of shared read counts. In case of purely unique mappings, taxa abundance estimates can be directly inferred from mapping counts.': ",sum(read.count.s2), " < ",N.s2, " (total number of reads in Sample 2)")   }
  }
  
  
  
 return (list("taxa"=taxa, "mat"=sim.mat, "counts_s1"= read.count.s1, "N1"=N.s1, "counts_s2"= read.count.s2, "N2"=N.s2) )

}   #end
