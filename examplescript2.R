params.groups <- newSplatParams(batchCells = 200, nGenes=1000) #parameters for control set
params.groups2 <- newSplatParams(batchCells = 200, nGenes=1000) #parameters for exp. cond. 1 set 
simdc1 <- splatSimulateGroups(params.groups, group.prob=c(0.25,0.25,0.5)) #simulation data for control set
simdc2 <- splatSimulateGroups(params.groups2, group.prob=c(0.15,0.15,0.7)) #simulation for the exp set 
simdc1norm <- normalize(simdc1)
plotPCA(simdc1norm)

#Count Data Creation 
countdc1 <- counts(simdc1) 
