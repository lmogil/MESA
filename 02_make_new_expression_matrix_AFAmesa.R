####Lauren S. Mogil####

#make new expression matrix###
###based on gencode v.18##


args <- commandArgs(trailingOnly=T)
#args <- c('22','0.5','1e6','AFA','no')
"%&%" = function(a,b) paste(a,b,sep="")

###############################################
### Directories & Variables

#wheelerlab2:
my.dir <- "/home/wheelerlab2/MESA_2017-04-01/"
lm.dir<-"/home/lauren/"

cohort <- args[1]  


################################################
### Functions & Libraries

library(glmnet)
library(data.table)
library(dplyr)

################################################
### Input files
expfile <- lm.dir %&%"meqtl_sorted_ALL_MESA_Epi_GEX_data_sidno.txt"
genelocfile <- lm.dir %&% "MESA_GEX_ilmn_gencode.v18.annotation.txt"

##MESA population expression files
#expfile <- my.dir %&%"heathers_meqtl_raw_exp/meqtl_sorted_AFA_MESA_Epi_GEX_data_sidno.txt"

#expfile <- my.dir %&%"heathers_meqtl_raw_exp/meqtl_sorted_CAU_MESA_Epi_GEX_data_sidno.txt"

#expfile <- my.dir %&%"heathers_meqtl_raw_exp/meqtl_sorted_HIS_MESA_Epi_GEX_data_sidno.txt"

################################################
### Read in expression data
expdata <- fread(expfile,header=T)
expmat <- data.matrix(expdata[,-1])
rownames(expmat) <- expdata$PROBE_ID
t.expmat <- t(expmat) #want genes in cols

##read gene location data with ILMN id and ENSG ids
geneloc <- data.frame(fread(genelocfile)) 
rownames(geneloc) <- geneloc[,1]
t.expmat <- t.expmat[,intersect(colnames(t.expmat),rownames(geneloc))] ###pull gene expression data w/gene info
                
expsamplelist <- rownames(t.expmat) ##samples with exp data##

reffilt<-geneloc[match(colnames(t.expmat),geneloc$PROBE_ID),] #match column names to ILMN ids

reffilt$gene_id<-as.character(reffilt$gene_id) 

colnames(t.expmat)[match(colnames(t.expmat),reffilt$PROBE_ID)]<- reffilt$gene_id #replace ILMN ids with ENSG ids


n_occur <- data.frame(table(colnames(t.expmat)))
n_occur_freq<-n_occur[n_occur$Freq > 1,] #find ENSG ids that occur more than 1 time
t.expmat_ens1<-t.expmat[,colnames(t.expmat) %in% n_occur$Var1[n_occur$Freq > 1]] #matrix of genes that occur more than once
t.expmat_ens2<-t.expmat[,colnames(t.expmat) %in% n_occur$Var1[n_occur$Freq == 1]] #matrix of genes that occur one time

##define row and column names
colnames = colnames(t.expmat_ens1)
rownames=rownames(t.expmat_ens1)

# calculate the means of the column names that match from the matrix that has genes that occur more than once

means <- data.frame(sapply(unique(colnames), function(name) rowMeans(t.expmat_ens1[,colnames(t.expmat_ens1) == name])))
#1163 2717+7426

##merge means with genes that occur once 
final_exp_mat<-cbind2(t.expmat_ens2,means)
#final number of genes 10143

t.final<-t(final_exp_mat)

write.table(t.final, '/home/lauren/ALL_MESA_EPI_GEX_ens.txt', quote = FALSE, sep = "\t")