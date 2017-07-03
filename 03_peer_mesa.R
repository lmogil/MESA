####Lauren S. Mogil###
###run PEER factors 10-100 for all pops 
"%&%" = function(a,b) paste(a,b,sep="")
library(peer)
library(tibble)
library(dplyr)
#loop through pops
for(pop in c('AFA','HIS','CAU','ALL')){
expr = read.table('/home/lauren/'%&% pop %&%'_MESA_EPI_GEX_ens.txt',header=T,check.names = F)
expmat = t(as.matrix(expr[,-1])) #need genes in cols, ids in rows
  ## i is number of peer factors to calculate, recommend 25% of sample size, but no more than 100
for(i in c(10,20,30,50,100)){
    model = PEER()
    PEER_setPhenoMean(model,expmat)
    PEER_setNk(model,i)
    PEER_update(model)
    factors = PEER_getX(model)
    weights = PEER_getW(model)
    precision = PEER_getAlpha(model)
    residuals = PEER_getResiduals(model)
    print(dim(residuals))
    cat(pop, "Nk =", i)
    print(PEER_plotModel(model))
    adjexp <- t(residuals) 
    print(dim(adjexp))
    rownames(adjexp) <- expr[,1]
    colnames(adjexp) <- colnames(expr)[-1]
    adjexp <- as.data.frame(adjexp) %>% rownames_to_column('PROBE_ID')
    png(filename = '/home/lauren/new_peer_ensexp/' %&% pop %&% '_MESA_Epi_GEX_data_sidno_Nk-' %&% i %&% '.png')
    PEER_plotModel(model)
    dev.off()
    write.table(adjexp, file='/home/lauren/new_peer_ensexp/' %&% pop %&% '_MESA_Epi_GEX_data_sidno_Nk-' %&% i %&% '.txt', quote=F, row.names=F, sep='\t')
}
}
