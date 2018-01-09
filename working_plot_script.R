

#########
PCA plot fig1 wheeler lab 2

c<-ggplot() + geom_point(data=gwas_cau,aes(x=PC1,y=PC2,col=pop,shape=pop))+geom_point(data=hm3_cau,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw()+ theme(text = element_text(size=20),legend.text=element_text(size=20)) + scale_colour_brewer(palette="Set1")
a<-ggplot() + geom_point(data=gwas_hw,aes(x=PC1,y=PC2,col=pop,shape=pop))+geom_point(data=hm3_hw,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw()+ theme(text = element_text(size=20),legend.text=element_text(size=20)) + scale_colour_brewer(palette="Set1")
b<-ggplot() + geom_point(data=gwas_h,aes(x=PC1,y=PC2,col=pop,shape=pop))+geom_point(data=hm3_h,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw() + theme(text = element_text(size=20),legend.text=element_text(size=20))+ scale_colour_brewer(palette="Set1")
d<-ggplot() + geom_point(data=gwas_c3,aes(x=PC1,y=PC2,col=pop,shape=pop))+ geom_point(data=gwas_y3,aes(x=PC1,y=PC2,col=pop,shape=pop))+ geom_point(data=gwas_asn3,aes(x=PC1,y=PC2,col=pop,shape=pop))+ geom_point(data=gwas_p3wp,aes(x=PC1,y=PC2,col=pop,shape=pop))+ theme_bw()+theme(text = element_text(size=20),legend.text=element_text(size=20)) + scale_colour_brewer(palette="Set1")

  

pdf("/home/lauren/Fig-mesa_peer.pdf",width=10,height=10)
multiplot(a+ggtitle('A\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),b + ggtitle('B\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),c + ggtitle('C\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),d + ggtitle('D\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols=2)
dev.off()


#########
FDR and matrix fig 2
for(i in 1:22){
  AFAres<-fread("zcat /home/lauren/MatrixEQTL_results/AFA_Nk_20_PFs_chr" %&% i %&% ".meqtl.cis.2017-07-05.txt.gz",header=T)
  CAUres<-fread("zcat /home/lauren/MatrixEQTL_results/CAU_Nk_20_PFs_chr" %&% i %&% ".meqtl.cis.2017-07-05.txt.gz",header=T)
  HISres<-fread("zcat /home/lauren/MatrixEQTL_results/HIS_Nk_20_PFs_chr" %&% i %&% ".meqtl.cis.2017-07-05.txt.gz",header=T)
  if(exists("AFAall")){
    AFAall <- rbind(AFAall,AFAres)
    CAUall <- rbind(CAUall,CAUres)
    HISall <- rbind(HISall,HISres)
  }else{
    AFAall <- AFAres
    CAUall <- CAUres
    HISall <- HISres
  }
}

#use each population as the discovery cohort, pull top hits (FDR < 0.05) and calc pi_1 in each pop as the replication cohort


poplist <- c('AFA','CAU','HIS')
npops <- length(poplist)
pi1_matrix <- matrix(NA,nrow=npops,ncol=npops)
rownames(pi1_matrix) <- poplist
colnames(pi1_matrix) <- poplist

for(i in 1:length(poplist)){
  for(j in 1:length(poplist)){
    if(poplist[i] == poplist[j]){
      pi1 <- 1
    }else{
      pop1 <- get(poplist[i] %&% "all")
      pop2 <- get(poplist[j] %&% "all")
      pop1fdr05 <- dplyr::filter(pop1, FDR < 0.05)
      pop2tested <- inner_join(pop1fdr05,pop2,by=c("snps","gene"))
      pop2pvals <- pop2tested$pvalue.y
      qobjpop2 <- qvalue(p = pop2pvals)
      print("Disc: " %&% poplist[i] %&% " Rep: " %&% poplist[j])
      print(hist(qobjpop2))
      pi1 <- 1 - qobjpop2$pi0
    }
    pi1_matrix[i,j] <- pi1
  }
}

print(pi1_matrix)

pi1_melt <- melt(pi1_matrix)
pi1_melt$Var1 = with(pi1_melt, factor(Var1, levels = rev(levels(Var1))))
s<-ggplot(pi1_melt, aes(x=Var2,y=Var1,fill=value)) + geom_raster() + xlab("Replication Population") + 
  ylab("Discovery Population") + theme_bw(25) + scale_fill_gradient(name=expression(pi[1])) 

mesasumahc<-read.table("/home/lauren/MatrixEQTL_results/MESA_meqtl.cis_summaryahc.txt", header = T)
mesasumahc$Nk <- as.numeric(mesasumahc$Nk)
l<-ggplot(mesasumahc,aes(x=Nk,y=FDR_0.05,col=pop,shape=pop)) + geom_point()+ geom_line()+theme_bw()+theme(text = element_text(size=15),legend.text=element_text(size=15)) + scale_colour_brewer(palette="Set1")+coord_cartesian(xlim=c(0,100),ylim=c(3.0e+5,2.50e+6))+ylab(expression("Number of eQTLs with FDR < 0.05"))+xlab('Number of PEER factors')
#vp <- venn.diagram(list(AFA=SNP_pop_A$SNP,HIS=SNP_pop_H$SNP,CAU=SNP_pop_C$SNP),fill = 2:4, alpha = 0.3, filename = NULL)
pdf("/home/lauren/MESA_fdr_mat.pdf",width=10,height=5)
multiplot(l+ggtitle('A\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),s+ggtitle('B\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols = 2)
dev.off()


##############
popcorn fig3

p<-ggplot(plot_this, aes(x = h2, y = Variable, color = Variable)) +
  geom_smooth(aes(x = h2, y = AFA_CAU, color = "AFA-CAU")) +
  geom_smooth(aes(x = h2, y = AFA_HIS, color = "AFA-HIS")) +
  geom_smooth(aes(x = h2, y = CAU_HIS, color = "CAU-HIS")) +
  scale_color_manual(name = "pop pair", values = c("blue", "magenta", "green")) + 
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "h2 threshold in AFA", y = "mean rG")
q<-ggplot(allrG_allh2)+ theme(plot.title = element_text(hjust = 0.5)) +geom_density(aes(x=rG_AFA_CAU, color = "AFA-CAU",fill="AFA-CAU"),color="blue" ,size=0.5, alpha = 1) + scale_fill_manual( values = c("blue","magenta","green"))+geom_density(aes(x=rG_AFA_HIS, color = "AFA-HIS",fill="AFA-HIS"),color="magenta" , size=0.5, alpha = 1) +geom_density(aes(x=rG_CAU_HIS, color = "CAU-HIS",fill="CAU-HIS"), color="green" ,size=0.5, alpha = 0.5) + geom_vline(xintercept = mean(allrG_allh2$rG_AFA_CAU, na.rm = T), color = "blue") +geom_vline(xintercept = mean(allrG_allh2$rG_AFA_HIS, na.rm = T), color = "magenta") + geom_vline(xintercept = mean(allrG_allh2$rG_CAU_HIS, na.rm = T), color = "green") + scale_color_manual(name = "MESA crosses", values = c("blue", "magenta", "green")) + labs(x = "rG")  
q<-q+labs(fill = "pop pair")
pdf("/home/lauren/Fig-gctacomp.pdf",width=10,height=5)
multiplot(q+ggtitle('A\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),p + ggtitle('B\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols=2)
dev.off()



####Corr plot MESA only fig5a
m<-ggplot(allR2,aes(AFA,HIS)) + geom_point()+theme_bw()+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+xlab(expression("AFA R"^2))+ylab(expression("HIS R"^2))
n<-ggplot(allR2,aes(CAU,AFA)) + geom_point()+theme_bw()+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+xlab(expression("CAU R"^2))+ylab(expression("AFA R"^2))
o<-ggplot(allR2,aes(HIS,CAU)) + geom_point()+theme_bw()+coord_cartesian(xlim=c(0,1),ylim=c(0,1))+xlab(expression("HIS R"^2))+ylab(expression("CAU R"^2))
pdf("/home/lauren/Fig-mesa-corr.pdf",width=8,height=3)
multiplot(m+ggtitle('A\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),n + ggtitle('B\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),o+ ggtitle('C\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols=3)
dev.off()


##### non-euro vs euro fig5b
g<-ggplot(allR2,aes(R2aR2c,AFA,colour=CAU)) + geom_point()+theme_bw()+coord_cartesian(xlim=c(-0.75,0.75),ylim=c(0,1))+ylab(expression("AFA R"^2))+xlab(expression("AFA R"^{2}~"-CAU R"^{2}))+scale_color_gradient(name= expression("CAU R"^2),low="blue", high="pink")+geom_vline(xintercept = 0.2)+geom_vline(xintercept = -0.2)

h<-ggplot(allR2,aes(R2hR2c,HIS,colour=CAU)) + geom_point()+theme_bw()+coord_cartesian(xlim=c(-0.75,0.75),ylim=c(0,1))+ylab(expression("HIS R"^2))+xlab(expression("HIS R"^{2}~"-CAU R"^{2}))+scale_color_gradient(name= expression("CAU R"^2),low="blue", high="pink")+geom_vline(xintercept = 0.2)+geom_vline(xintercept = -0.2)

pdf("/home/lauren/Fig-MESA-COMP-POP.pdf",width=8,height=3)
multiplot(g+ggtitle('D\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),h + ggtitle('E\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols=2)
dev.off()

pdf("/home/lauren/Fig-mesa-corr_comp.pdf",width=8,height=3)
multiplot(m+ggtitle('A\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),n + ggtitle('B\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),o+ ggtitle('C\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18))+(g+ggtitle('D\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),h + ggtitle('E\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols=3)
dev.off()




##############R2lasso v alpha fig 4
alphah1 <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/HIS_all_Results_1.txt',header=TRUE) %>%  mutate(`1`=in_sample_R2) %>% select(gene_id,`1`)
ngenesallh <- length(unique(alpha1$gene))
ngenesallh
alpha50h <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/old_results/HIS_all_Results.txt',header=TRUE) %>% mutate(`0.50`=in_sample_R2) %>% select(gene_id,`0.50`)
alpha0h <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/HIS_all_Results_0.05.txt',header=TRUE) %>% mutate(`0.05`= in_sample_R2) %>% select(gene_id,`0.05`)

datah <- inner_join(alpha0h,alpha50h,by='gene_id')
datah<- inner_join(datah,alphah1,by='gene_id')
gdatah <- gather(datah,alpha,R2,2:3)
gdatah$pop<-"HIS"

alphaa1 <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/AFA_all_Results_1.txt',header=TRUE) %>%  mutate(`1`=in_sample_R2) %>% select(gene_id,`1`)
ngenesalla <- length(unique(alpha1$gene))
ngenesalla
alpha50a <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/old_results/AFA_all_Results.txt',header=TRUE) %>% mutate(`0.50`=in_sample_R2) %>% select(gene_id,`0.50`)
alpha0a <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/AFA_all_Results_0.05.txt',header=TRUE) %>% mutate(`0.05`= in_sample_R2) %>% select(gene_id,`0.05`)

dataa <- inner_join(alpha0a,alpha50a,by='gene_id')
dataa<- inner_join(dataa,alphaa1,by='gene_id')
gdataa <- gather(dataa,alpha,R2,2:3)
gdataa$pop<-"AFA"

alpha1 <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/CAU_all_Results_1.txt',header=TRUE) %>%  mutate(`1`=in_sample_R2) %>% select(gene_id,`1`)
ngenesall <- length(unique(alpha1$gene))
ngenesall
alpha50 <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/old_results/CAU_all_Results.txt',header=TRUE) %>% mutate(`0.50`=in_sample_R2) %>% select(gene_id,`0.50`)
alpha0 <- read.table('~/PredictDB_Pipeline_GTEx_v7/new_output/CAU_all_Results_0.05.txt',header=TRUE) %>% mutate(`0.05`= in_sample_R2) %>% select(gene_id,`0.05`)

data <- inner_join(alpha0,alpha50,by='gene_id')
data<- inner_join(data,alpha1,by='gene_id')
gdata <- gather(data,alpha,R2,2:3)
gdata$pop<-"CAU"


gdataahc<-full_join(gdata,gdatah)
gdataahc<-full_join(gdataahc,gdataa)

q<-ggplot(gdataahc, aes(y = `1` - R2, x = `1`, group=alpha, color=alpha)) + geom_point(show.legend = TRUE) + ylab(expression(paste("R"^2, " difference (LASSO - alpha)"))) + xlab(expression(paste("R"^2, " (LASSO)"))) +theme_bw(12)+ theme(legend.justification=c(0.05,1), legend.position=c(0,1))+coord_cartesian(xlim=c(0,1),ylim=c(-1,1))
pdf("/home/lauren/Fig-en_mesa.pdf",width=8,height=3)
q + facet_wrap( ~ pop, ncol=3)
dev.off()

##########hap r vs r2 fig 6


y<-ggplot(RvsR2ym_mesa) + geom_smooth(aes(x = ALL, y = MEX_all, color = "ALL")) +geom_smooth(aes(x = CAU, y = MEX_cau, color = "CAU")) +geom_smooth(aes(x = HIS, y = MEX_his, color = "HIS")) +geom_smooth(aes(x = AFA, y = MEX_afa, color = "AFA")) +scale_color_manual(name = "Model", values = c("blue", "magenta", "green","purple","red")) + theme(plot.title = element_text(hjust = 0.5)) +labs(x = "Model R2", y = "MXL R")

z<-ggplot(RvsR2ym_mesa) + geom_smooth(aes(x = ALL, y = YRI_all, color = "ALL")) +geom_smooth(aes(x = CAU, y = YRI_cau, color = "CAU")) +geom_smooth(aes(x = HIS, y = YRI_his, color = "HIS")) +geom_smooth(aes(x = AFA, y = YRI_afa, color = "AFA")) +scale_color_manual(name = "Model", values = c("blue", "magenta", "green","purple","red")) + theme(plot.title = element_text(hjust = 0.5)) +labs(x = "Model R2", y = "YRI R")


pdf("/home/lauren/Fig-rvr2.pdf",width=8,height=3)
multiplot(y+ggtitle('A\n')+theme(plot.title=element_text(hjust=0),text=element_text(size=18)),z + ggtitle('B\n')+ theme(plot.title=element_text(hjust=0),text=element_text(size=18)),cols=2)
dev.off()
