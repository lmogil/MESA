library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
"%&%" = function(a,b) paste(a,b,sep="")

px.dir = "/home/lauren/" 
obs.dir = "/home/lauren/"


pops <- c('MEX','YRI")
dbs <- c('AFA','HIS','CAU','AFHI','ALL')


for(d in dbs){
  for(pop in pops){
    
    predexp1 <- data.frame(fread(px.dir %&% "new_predixcan_mesa/"%&% pop %&%"_"%&% d %&%"_predicted_expression.txt"))
    rownames(predexp1) <- predexp1[,1]
    obsexp <- data.frame(fread(obs.dir %&% "geuvadis_data/CEU_expression_geu_Nk-10.txt"))
    #obsexp <- data.frame(fread(obs.dir %&% "HapMap_data/"%&% pop %&%"_exp_nk10_ens.txt"))
    rownames(obsexp)<-obsexp[,1]
    tobsexp <- t(obsexp[,-1]) #transpose the observed exp matrix
    
    #get the same genes in obs & pred and sort by ID and gene
    #obs2 <- data.frame(tobsexp[,colnames(tobsexp) %in% colnames(predexp1)])
    #obs1<- data.frame(obs2[,colnames(obs2) %in% afaalphafilt$gene_id])
    #obs <- obs1[order(rownames(obs1)),order(colnames(obs1))]
    
    obs2 <- data.frame(tobsexp[,colnames(tobsexp) %in% colnames(predexp1)])
    obs <- obs2[order(rownames(obs2)),order(colnames(obs2))]
    
    
    
    pred2 <- predexp1[,colnames(predexp1) %in% colnames(obs)]
    pred <- pred2[order(rownames(pred2)),order(colnames(pred2))]
    #convert to matrix and transpose
    predexp <- as.matrix(pred)
    obsexp <- as.matrix(obs)

    popres <- matrix(NA,ncol=1,nrow=dim(obsexp)[2])

    for(i in 1:dim(obsexp)[2]){
      corres <- cor.test(predexp[,i] , obsexp[,i])
      r <- signif(corres$estimate,3)
      popres[i,] <- r
    }
    if(exists("allres") == FALSE){
      allres = popres
    }else{
      allres<- cbind(allres,popres)
    }
  }

  colnames(allres) <- pops
  print(summary(allres))
  gres <- gather(data.frame(allres),key=pop,value=R)
  ggplot(gres,aes(x=pop,y=R,color=pop)) + geom_boxplot(outlier.size = 0.5) + theme_bw(15) + guides(color=FALSE) + ggtitle("Weights: " %&% d) + xlab("HapMap3 population")+ylab("Pearson's R (pred. v. obs.)")
  ggsave("/home/lauren/pvo_newexp/R_pred_v_obs_" %&% d %&% "_"%&% pop %&%"_db.png")
  rownames(allres) <- colnames(obs)
  write.table(allres,"/home/lauren/pvo_newexp/R_pred_v_obs_" %&% d %&% "_" %&% pop %&% "_db.txt",quote=F)
  rm("allres")
}

