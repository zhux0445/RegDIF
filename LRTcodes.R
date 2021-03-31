library(MASS)
library(mirt) 
mirtCluster(4)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
setwd('/Users/zhux0445/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("RESP1.csv",row.names = 1)
beta.1=read.csv("Beta1_ 1.dms",row.names = 1)
ad.1=read.csv("ADmat1_ 1.dms",row.names = 1)

soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }
# Dataset #4 (2 Non-uniform DIF items per scale)
J=20

N1=N2=N3=500 
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
N=N1+N2+N3

m=2
r=2
y1=c(0,0)
y2=c(1,0)
y3=c(0,1)

Amat1=params[,1:2]
Amat2=params[,3:4]
Amat3=params[,5:6]
Dmat1=matrix(params[,7],J,(m-1))
Dmat2=matrix(params[,8],J,(m-1))
Dmat3=matrix(params[,9],J,(m-1))
Amat1=as.matrix(Amat1)
Amat2=as.matrix(Amat2)
Amat3=as.matrix(Amat3)
rownames(Amat1) <- c()
rownames(Amat2) <- c()
rownames(Amat3) <- c()
#test cor=0.25
rep=1
set.seed(rep*100)
Theta=mvrnorm(n=N,mu=c(0,0),Sigma = matrix(c(1,0.25,0.25,1),2,2))
Theta1=Theta[1:N1,]
Theta2=Theta[(N1+1):(N1+N2),]
Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
resp=rbind(datasetR,datasetF1,datasetF2)


#####  No anchor


mirt.p.mat1=matrix(0,50,20) # 16 is the number of replications
power1=numeric(50)
tpI1=numeric(50)
bias.mirt1=matrix(0,50,3)
rmse.mirt1=matrix(0,50,3)
difrec.mirt.fn1=matrix(0,50,3) # DIF magnitude recovery with false negative included (only need the first column)
difrec.mirt.fn1.r=matrix(0,50,3) # DIF magnitude recovery for regularization method omnibus (only need the first column)
mirt.p.mat2=matrix(0,50,20) # 16 is the number of replications
power2=numeric(50)
tpI2=numeric(50)
bias.mirt2=matrix(0,50,3)
rmse.mirt2=matrix(0,50,3)
difrec.mirt.fn2=matrix(0,50,3) #(only need the second column)
mirt.p.mat3=matrix(0,50,20) # 16 is the number of replications
power3=numeric(50)
tpI3=numeric(50)
bias.mirt3=matrix(0,50,3)
rmse.mirt3=matrix(0,50,3)
difrec.mirt.fn3=matrix(0,50,3) #(only need the third column)

for (rep in 2:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('slopes'))
  dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add_sequential',items2test=c(1:J))
  power1[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif1))
  tpI1[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif1))
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('slopes',colnames(resp)[which(colnames(resp)%in%rownames(dif1)==0)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-c(4,5,12,13)]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[rep,2:3]= c(mean(abs((coef(md.refit0,simplify=T)$G2$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit0,simplify=T)$G3$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3])) #omnibus DIF recovery in report (only need this)
  #difrec.mirt.fn1.r[rep,2:3]= c(mean(abs((coef(md.refit.r,simplify=T)$G2$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit.r,simplify=T)$G3$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  #difrec.mirt.fn1.r[rep,1]=mean(c(difrec.mirt.fn1.r[rep,2],difrec.mirt.fn1.r[rep,3])) #omnibus DIF recovery in report (only need this)
  
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('slopes'))
  dif2=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add_sequential',items2test=c(1:J))
  power2[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif2))
  tpI2[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif2))
  if ((power2[rep])==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('slopes',colnames(resp01)[which(colnames(resp01)%in%rownames(dif2)==0)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('slopes'))
  dif3=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add_sequential',items2test=c(1:J))
  power3[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif3))
  tpI3[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif3))
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('slopes',colnames(resp02)[which(colnames(resp02)%in%rownames(dif3)==0)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1))
}
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)
