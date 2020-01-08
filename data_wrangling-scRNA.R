#1Loading required packages 


library(Seurat)
library(dplyr)
library(Matrix)
library(gdata)

#Comment: successfully installed all the libraries

#2 Reading NSCLC counts matrix 
#2a Identifying the pathname: /Users/hephaeschuen/Documents/GitHub/methods-integrative-analysis-scrnaseq-bbms4001/filtered_gene_bc_matrices/GRCh38/matrix.mtx
counts <- Read10X(data.dir = "/Users/hephaeschuen/Documents/GitHub/methods-integrative-analysis-scrnaseq-bbms4001/filtered_gene_bc_matrices/GRCh38")
#problem: What is dgCMatrix class? What are its properties?
#problem: Is it possible to obtain the count matrix in the dgCMatrix, and convert it
#to standard 2D matrix?
#answer: Yes it is possible. Use as.matrix(). But the size of the matrix in the memory 
#will be extremely large 

#3 Summary Statistics 
counts_per_cell <- Matrix::colSums(counts)
counts_per_gene <- Matrix::rowSums(counts)
genes_per_cell <- Matrix::colSums(counts>0) # count gene only if it has non-zero reads mapped.



