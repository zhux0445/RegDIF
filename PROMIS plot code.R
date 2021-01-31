library(ggplot2)
item=seq(1,21,1)
uniformLRT1=c(0,.481,0,-.071,0,.159,.282,0,0,0,-.455,-.388,-.347,0,0,0,0,-.423,0,0,0)
uniformLRT2=c(0,.729,0,-.516,0,-.258,-.079,0,0,0,-.774,-.507,-.722,0,0,0,0,-.646,0,0,0)
uniformEM1=c(rep(0,10),-.351,rep(0,10))
uniformEM2=c(rep(0,3),-.39,rep(0,6),-.672,0,-.463,rep(0,4),-.339,rep(0,3))
uniformEMM1=c(rep(0,10),-.35,rep(0,10))
uniformEMM2=c(rep(0,3),-.55,0,-.416,-.328,rep(0,3),-.672,0,-.465,rep(0,4),-.339,rep(0,3))
uniformAdapt1=c(rep(0,3),-.24,rep(0,6),-.351,rep(0,10))
uniformAdapt2=c(rep(0,3),-.606,0,-.347,rep(0,4),-.673,0,-.463,rep(0,4),-.340,rep(0,3))
uniform=cbind(item,uniformLRT1,uniformLRT2,uniformEM1,uniformEM2,uniformEMM1,uniformEMM2,uniformAdapt1,uniformAdapt2)
uniform=as.data.frame(uniform)
p <- ggplot(uniform, aes(item,uniformLRT1))
p + theme_bw()+  geom_point(uniform,  mapping=aes(x=as.factor(item),y=uniformLRT1),colour="gray",shape=15,size=3) +  geom_point(uniform,  mapping=aes(x=item,y=uniformLRT2),shape=15,size=3) +  geom_point(uniform,  mapping=aes(x=item,y=uniformEM1),colour="gray",shape=16,size=3) +  geom_point(uniform,  mapping=aes(x=item,y=uniformEM2),shape=16,size=3) +
  geom_point(uniform,  mapping=aes(x=item,y=uniformEMM1),colour="gray",shape=17,size=3) +  geom_point(uniform,  mapping=aes(x=item,y=uniformEMM2),shape=17,size=3) +
  geom_point(uniform,  mapping=aes(x=item,y=uniformAdapt1),colour="gray",shape=18,size=3) +  geom_point(uniform,  mapping=aes(x=item,y=uniformAdapt2),shape=18,size=3) 



uniform2=cbind(rep(item,8),c(uniformLRT1,uniformLRT2,uniformEM1,uniformEM2,uniformEMM1,uniformEMM2,uniformAdapt1,uniformAdapt2),c(rep("LRT",42),rep("Lasso EM",42),rep("Lasso EMM",42),rep( "Adaptive Lasso",42)),rep(c(rep("Age 50-64",21),rep("Age 65-84",21)),4))
uniform2=as.data.frame(uniform2)
colnames(uniform2)=c("Item","DIF_magnitude","Method","Group")

uniform2$Item=factor(uniform2$Item,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))
uniform2$DIF_magnitude[uniform2$DIF_magnitude==0]=NA
p <- ggplot(na.omit(uniform2), aes(x=Item,y=as.numeric(DIF_magnitude),shape=Method,color=Group))
p+geom_point(size=4)+ scale_color_grey()+theme_bw(base_size=20)+xlab("Item")+ylab("DIF Magnitude Estimates")+ geom_abline(intercept = 0, slope = 0)

NonuniformLRT1=rep(0,21)
NonuniformLRT2=rep(0,21)
NonuniformEM1=c(rep(0,10),.653,.658,.530,rep(0,4),.94,rep(0,3))
NonuniformEM2=c(rep(0,3),.572,0,.632,rep(0,4),1.295,.736,.951,rep(0,4),1.024,rep(0,3))
NonuniformEMM1=c(rep(0,10),.648,.661,.528,rep(0,4),.945,rep(0,3))
NonuniformEMM2=c(rep(0,3),.797,0,.856,.612,rep(0,3),1.297,.742,.958,rep(0,4),1.021,rep(0,3))
NonuniformAdapt1=c(rep(0,10),.650,.659,.528,rep(0,4),.942,rep(0,3))
NonuniformAdapt2=c(rep(0,3),.785,0,.843,.601,rep(0,3),1.301,.745,.962,rep(0,4),1.028,rep(0,3))
Nonuniform=cbind(rep(item,8),c(NonuniformLRT1,NonuniformLRT2,NonuniformEM1,NonuniformEM2,NonuniformEMM1,NonuniformEMM2,NonuniformAdapt1,NonuniformAdapt2),c(rep("LRT",42),rep("Lasso EM",42),rep("Lasso EMM",42),rep( "Adaptive Lasso",42)),rep(c(rep("Age 50-64",21),rep("Age 65-84",21)),4))
Nonuniform=as.data.frame(Nonuniform)
colnames(Nonuniform)=c("Item","DIF_magnitude","Method","Group")

Nonuniform$Item=factor(Nonuniform$Item,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))
Nonuniform$DIF_magnitude[Nonuniform$DIF_magnitude==0]=NA
p2 <- ggplot(na.omit(Nonuniform), aes(x=Item,y=as.numeric(DIF_magnitude),shape=Method,color=Group))
p2+geom_point(size=4)+ scale_color_grey()+theme_bw(base_size=20)+xlab("Item")+ylab("DIF Magnitude Estimates")


OmniNonuniformLRT1=c(0,-.419,rep(0,8),.077,.18,.079,0,0,-.591,rep(0,5))
OmniNonuniformLRT2=c(0,-.375,rep(0,8),.469,.016,-.073,0,0,-.861,rep(0,5))
OmniNonuniformEM1=rep(0,21)
OmniNonuniformEM2=rep(0,21)
OmniNonuniformEMM1=rep(0,21)
OmniNonuniformEMM2=rep(0,21)
OmniNonuniformAdapt1=rep(0,21)
OmniNonuniformAdapt2=c(rep(0,15),1.27,rep(0,5))
OmniUniformLRT1=c(0,.505,rep(0,8),-.373,-.244,-.209,0,0,.111,rep(0,5))
OmniUniformLRT2=c(0,.857,rep(0,8),-.612,-.416,-.705,0,0,.578,rep(0,5))
OmniUniformEM1=c(rep(0,10),-.35,rep(0,10))
OmniUniformEM2=c(.391,.609,rep(0,8),-.646,0,-.435,rep(0,8))
OmniUniformEMM1=c(0,.492,rep(0,8),-.35,rep(0,10))
OmniUniformEMM2=c(.391,.915,rep(0,8),-.672,0,-.462,rep(0,4),-.34,rep(0,3))
OmniUniformAdapt1=c(0,.493,rep(0,19))
OmniUniformAdapt2=c(.368,.889,0,-.313,0,.843,.601,rep(0,3),-.417,0,-.433,rep(0,2),.849,rep(0,5))
Omni=cbind(rep(item,16),c(OmniNonuniformLRT1,OmniNonuniformLRT2,OmniNonuniformEM1,OmniNonuniformEM2,OmniNonuniformEMM1,OmniNonuniformEMM2,OmniNonuniformAdapt1,OmniNonuniformAdapt2,OmniUniformLRT1,OmniUniformLRT2,OmniUniformEM1,OmniUniformEM2,OmniUniformEMM1,OmniUniformEMM2,OmniUniformAdapt1,OmniUniformAdapt2),rep(c(rep("LRT",42),rep("Lasso EM",42),rep("Lasso EMM",42),rep( "Adaptive Lasso",42)),2),rep(c(rep("Age 50-64",21),rep("Age 65-84",21)),8),c(rep("DIF on Slope",168),rep("DIF on Intercept",168)))
Omni=as.data.frame(Omni)
colnames(Omni)=c("Item","DIF_magnitude","Method","Group","DIF_type")

Omni$Item=factor(Omni$Item,levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21"))
Omni$DIF_magnitude[Omni$DIF_magnitude==0]=NA
p3 <- ggplot(na.omit(Omni), aes(x=Item,y=as.numeric(DIF_magnitude),shape=Method,color=DIF_type))
p3+geom_point(size=4)+ scale_color_grey()+theme_bw(base_size=20)+xlab("Item")+ylab("DIF Magnitude Estimates")+facet_grid(.~Group)+ geom_abline(intercept = 0, slope = 0)


