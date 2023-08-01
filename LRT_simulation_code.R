library(MASS)
library(mirt) 
mirtCluster(8)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
setwd('/Users/zhux0445/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("RESP1new_SS600.csv",row.names = 1)
beta.1=read.csv("Beta1_ 1.dms",row.names = 1)
ad.1=read.csv("ADmat1_ 1.dms",row.names = 1)

soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }
# Dataset #4 (2 Non-uniform DIF items per scale)
J=20

N1=N2=N3=200 
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


#####  No anchor (program Carol Woods select procedure)


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
anchors.all=matrix(0,50,2)
for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #select anchor
  md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
  d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop')
  ratio1=d$X2/d$df
  anchor1=which(ratio1==sort((d$X2/d$df)[c(1,3:11)])[1])
  anchor2=which(ratio1==sort((d$X2/d$df)[c(2,12:20)])[1])
  anchors.all[rep,]=c(anchor1,anchor2)
  print(anchors.all[rep,])
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor1,anchor2)]))#,anchor3,anchor4)]))
  dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor1,anchor2)])#,anchor3,anchor4)])
  dif1.t=dif1[which(dif1$adj_pvals<0.05),]
  power1[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif1.t))
  tpI1[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif1.t))
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
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
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor1,anchor2)]))#,anchor3,anchor4)]))
  dif2=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor1,anchor2)])#,anchor3,anchor4)])
  dif2.t=dif2[which(dif2$adj_pvals<0.05),]
  power2[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif2.t))
  tpI2[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif2.t))
  if ((power2[rep])==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[which(colnames(resp01)%in%rownames(dif2.t)==0)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[c(anchor1,anchor2)]))#,anchor3,anchor4)]))
  dif3=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-c(anchor1,anchor2)])#,anchor3,anchor4)])
  dif3.t=dif3[which(dif3$adj_pvals<0.05),]
  power3[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif3.t))
  tpI3[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif3.t))
  if ((power3[rep])==0){
    md.refit02 <-multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:J]))
  } else {
    md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[which(colnames(resp02)%in%rownames(dif3.t)==0)]))
  }
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1))
}
sum(power1)/(50*12)
sum(tpI1)/(50*8)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(power2)/(50*4)
sum(tpI2)/(50*16)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

sum(power3)/(50*4)
sum(tpI3)/(50*16)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(difrec.mirt.fn3)



#####  No anchor (program iterative select procedure)


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
anchors.all=matrix(0,50,20)
for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #select anchor
  anchors0=c(1:20)
  diff=1
  while(diff>0){
    md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors0]))
    d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop')
    anchors=which(d$adj_pvals>0.05)
    diff=length(anchors0)-length(anchors)
    anchors0=anchors
    #anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
  }
  anchors.all[rep,1:length(anchors)]=anchors#,anchor3,anchor4)
  print(anchors.all[rep,1:length(anchors)])
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])#,anchor3,anchor4)])
  dif1.t=dif1[which(dif1$adj_pvals<0.05),]
  #power1[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif1.t))
  #tpI1[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif1.t))
  power1[rep]=sum(c("V4" , "V5" ,"V12", "V13")%in%rownames(dif1.t))
  tpI1[rep]=sum(c("V1","V2","V3", "V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif1.t))
  
  
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif2=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])#,anchor3,anchor4)])
  dif2.t=dif2[which(dif2$adj_pvals<0.05),]
  #power2[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif2.t))
  #tpI2[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif2.t))
  power2[rep]=sum(c("V4" , "V5" ,"V12", "V13")%in%rownames(dif2.t))
  tpI2[rep]=sum(c("V1","V2","V3", "V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif2.t))
  if ((power2[rep])==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[which(colnames(resp01)%in%rownames(dif2.t)==0)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  z3=(coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]
  difrec.mirt.fn2[rep,2]= mean(abs(z3[z3!=0]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif3=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])#,anchor3,anchor4)])
  dif3.t=dif3[which(dif3$adj_pvals<0.05),]
  #power3[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif3.t))
  #tpI3[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif3.t))
  power3[rep]=sum(c("V4" , "V5" ,"V12", "V13")%in%rownames(dif3.t))
  tpI3[rep]=sum(c("V1","V2","V3", "V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif3.t))
  if ((power3[rep])==0){
    md.refit02 <-multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:J]))
  } else {
    md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[which(colnames(resp02)%in%rownames(dif3.t)==0)]))
  }
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  z4=(coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]
  difrec.mirt.fn3[rep,3]=mean(abs(z4[z4!=0]-1))
}
sum(power1)/(50*4)
sum(tpI1)/(50*16)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(na.omit(difrec.mirt.fn1))

sum(power2)/(50*4)
sum(tpI2)/(50*16)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(na.omit(difrec.mirt.fn2))


sum(power3)/(50*4)
sum(tpI3)/(50*16)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(na.omit(difrec.mirt.fn3))

write.csv(power1, file = "NALRTUnifpower1lowcor.csv")
write.csv(tpI1, file = "NALRTUniftpI1lowcor.csv")
write.csv(bias.mirt1, file = "NALRTUnifbias1lowcor.csv")
write.csv(rmse.mirt1, file = "NALRTUnifrmse1lowcor.csv")
write.csv(difrec.mirt.fn1, file = "NALRTUnifrec1lowcor.csv")
write.csv(power2, file = "NALRTUnifpower12lowcor.csv")
write.csv(tpI2, file = "NALRTUniftpI12lowcor.csv")
write.csv(difrec.mirt.fn2, file = "NALRTUnifrec12lowcor.csv")
write.csv(power3, file = "NALRTUnifpower13lowcor.csv")
write.csv(tpI3, file = "NALRTUniftpI13lowcor.csv")
write.csv(difrec.mirt.fn3, file = "NALRTUnifrec13lowcor.csv")
write.csv(anchors.all, file = "NALRTUnifanchor1lowcor.csv")

sum(read.csv("NALRTUnifpower1.csv")[,2])/(50*4)
sum(read.csv("NALRTUniftpI1.csv")[,2])/(50*16)
sum(read.csv("NALRTUnifpower12.csv")[,2])/(50*4)
sum(read.csv("NALRTUniftpI12.csv")[,2])/(50*16)
sum(read.csv("NALRTUnifpower13.csv")[,2])/(50*4)
sum(read.csv("NALRTUniftpI13.csv")[,2])/(50*16)



sum(power1)/(50*12)
sum(tpI1)/(50*8)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(power2)/(50*12)
sum(tpI2)/(50*8)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

sum(power3)/(50*12)
sum(tpI3)/(50*8)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(difrec.mirt.fn3)


#####  No anchor (drop_sequential select)


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
anchors.all=matrix(0,50,20)
for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #select anchor
  md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
  d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop_sequential')
  anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
  anchors.all[rep,1:length(anchors)]=anchors#,anchor3,anchor4)
  print(anchors.all[rep,1:length(anchors)])
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])#,anchor3,anchor4)])
  dif1.t=dif1[which(dif1$adj_pvals<0.05),]
  power1[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif1.t))
  tpI1[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif1.t))
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
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
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif2=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])#,anchor3,anchor4)])
  dif2.t=dif2[which(dif2$adj_pvals<0.05),]
  power2[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif2.t))
  tpI2[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif2.t))
  if ((power2[rep])==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[which(colnames(resp01)%in%rownames(dif2.t)==0)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif3=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])#,anchor3,anchor4)])
  dif3.t=dif3[which(dif3$adj_pvals<0.05),]
  power3[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif3.t))
  tpI3[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif3.t))
  if ((power3[rep])==0){
    md.refit02 <-multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:J]))
  } else {
    md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[which(colnames(resp02)%in%rownames(dif3.t)==0)]))
  }
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1))
}
sum(power1)/(50*4)
sum(tpI1)/(50*16)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(power2)/(50*4)
sum(tpI2)/(50*16)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

sum(power3)/(50*4)
sum(tpI3)/(50*16)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(difrec.mirt.fn3)

#####  drop_sequential


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
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
  dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop_sequential')
  power1[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif1))
  tpI1[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif1))
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[which(colnames(resp)%in%rownames(dif1)==0)]))
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
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
  dif2=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop_sequential')
  power2[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif2))
  tpI2[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif2))
  if ((power2[rep])==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[which(colnames(resp01)%in%rownames(dif2)==0)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)))
  dif3=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop_sequential')
  power3[rep]=sum(c("V4" , "V5" , "V12", "V13")%in%rownames(dif3))
  tpI3[rep]=sum(c("V1","V2","V3","V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif3))
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[which(colnames(resp02)%in%rownames(dif3)==0)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1))
}
sum(power1)/(50*4)
sum(tpI1)/(50*16)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(power2)/(50*4)
sum(tpI2)/(50*16)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

sum(power3)/(50*4)
sum(tpI3)/(50*16)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(difrec.mirt.fn3)


sd(power1/4)/7
sd(tpI1/16)/7
sd(power2/4)/7
sd(tpI2/16)/7
sd(power3/4)/7
sd(tpI3/16)/7




#####  Non-uniform (dif on both slope and intercept)
#####  No anchor (program iterative select procedure)


mirt.p.mat1=matrix(0,50,20) # 16 is the number of replications
power1=numeric(50)
tpI1=numeric(50)
bias.mirt1=matrix(0,50,3)
rmse.mirt1=matrix(0,50,3)
difrec.mirt.fn1=matrix(0,50,3) # DIF magnitude recovery with false negative included (only need the first column)
difrec.mirt.fn1.a=matrix(0,50,3) 
#difrec.mirt.fn1.r=matrix(0,50,3) # DIF magnitude recovery for regularization method omnibus (only need the first column)
mirt.p.mat2=matrix(0,50,20) # 16 is the number of replications
power2=numeric(50)
tpI2=numeric(50)
bias.mirt2=matrix(0,50,3)
rmse.mirt2=matrix(0,50,3)
difrec.mirt.fn2=matrix(0,50,3) #(only need the second column)
difrec.mirt.fn2.a=matrix(0,50,3) 
mirt.p.mat3=matrix(0,50,20) # 16 is the number of replications
power3=numeric(50)
tpI3=numeric(50)
bias.mirt3=matrix(0,50,3)
rmse.mirt3=matrix(0,50,3)
difrec.mirt.fn3=matrix(0,50,3) #(only need the third column)
difrec.mirt.fn3.a=matrix(0,50,3)
anchors.all=matrix(0,50,20)
for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #select anchor
  anchors0=c(1:20)
  diff=1
  while(diff>0){
    md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors0]))
    d=DIF(md.cons0, which.par = c('a1','a2','d'), p.adjust = 'fdr',scheme = 'drop')
    anchors=which(d$adj_pvals>0.05)
    diff=length(anchors0)-length(anchors)
    anchors0=anchors
    #anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
  }
  anchors.all[rep,1:length(anchors)]=anchors#,anchor3,anchor4)
  print(anchors.all[rep,1:length(anchors)])
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif11=DIF(md.noncons0, which.par = c('a1','d'), p.adjust = 'fdr',scheme = 'add',items2test=(c(1:J)[-anchors])[which(c(1:J)[-anchors]%in%c(1,3:11))])
  dif12=DIF(md.noncons0, which.par = c('a2','d'), p.adjust = 'fdr',scheme = 'add',items2test=(c(1:J)[-anchors])[which(c(1:J)[-anchors]%in%c(2,12:20))])
  dif1.t=rbind(dif11[which(dif11$adj_pvals<0.05),],dif12[which(dif12$adj_pvals<0.05),])
  #power1[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif1.t))
  #tpI1[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif1.t))
  power1[rep]=sum(c("V4" , "V5" ,"V12", "V13")%in%rownames(dif1.t))
  tpI1[rep]=sum(c("V1","V2","V3", "V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif1.t))
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-c(4,5,12,13)]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  z1=(coef(md.refit0,simplify=T)$G2$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4:9,12:17)] # 4:9,12:17
  z2=(coef(md.refit0,simplify=T)$G3$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4:9,12:17)]
  difrec.mirt.fn1[rep,2:3]= c(mean(abs(z1[z1!=0]-0.25)),mean(abs(z2[z2!=0]-0.6)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3])) #omnibus DIF recovery in report (only need this)
  z3=(coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[c(4:9)]
  z4=(coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[c(12:17)]
  z5=(coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[c(4:9)]
  z6=(coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[c(12:17)]
  difrec.mirt.fn1.a[rep,2]= mean(c(abs(z3[z3!=0]-0.4),abs(z4[z4!=0]-0.4)))
  difrec.mirt.fn1.a[rep,3]= mean(c(abs(z5[z5!=0]-0.6),abs(z6[z6!=0]-0.6)))
  
  #difrec.mirt.fn1.r[rep,2:3]= c(mean(abs((coef(md.refit.r,simplify=T)$G2$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit.r,simplify=T)$G3$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  #difrec.mirt.fn1.r[rep,1]=mean(c(difrec.mirt.fn1.r[rep,2],difrec.mirt.fn1.r[rep,3])) #omnibus DIF recovery in report (only need this)
  
  
  
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif21=DIF(md.noncons01, which.par = c('a1','d'), p.adjust = 'fdr',scheme = 'add',items2test=(c(1:J)[-anchors])[which(c(1:J)[-anchors]%in%c(1,3:11))])#,anchor3,anchor4)])
  dif22=DIF(md.noncons01, which.par = c('a2','d'), p.adjust = 'fdr',scheme = 'add',items2test=(c(1:J)[-anchors])[which(c(1:J)[-anchors]%in%c(2,12:20))])#,anchor3,anchor4)])
  dif2.t=rbind(dif21[which(dif21$adj_pvals<0.05),],dif22[which(dif22$adj_pvals<0.05),])
  #power2[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif2.t))
  #tpI2[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif2.t))
  power2[rep]=sum(c("V4" , "V5" ,"V12", "V13")%in%rownames(dif2.t))
  tpI2[rep]=sum(c("V1","V2","V3", "V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif2.t))
  if ((power2[rep])==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp01)[which(colnames(resp01)%in%rownames(dif2.t)==0)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  z12=(coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4:9,12:17)] 
  difrec.mirt.fn2[rep,2]= mean(abs(z12[z12!=0]-0.25))
  z32=(coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[c(4:9)]
  z42=(coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[c(12:17)]
  difrec.mirt.fn2.a[rep,2]= mean(c(abs(z32[z32!=0]-0.4),abs(z42[z42!=0]-0.4)))
  
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors]))#,anchor3,anchor4)]))
  dif31=DIF(md.noncons02, which.par = c('a1','d'), p.adjust = 'fdr',scheme = 'add',items2test=(c(1:J)[-anchors])[which(c(1:J)[-anchors]%in%c(1,3:11))])#,anchor3,anchor4)])
  dif32=DIF(md.noncons02, which.par = c('a2','d'), p.adjust = 'fdr',scheme = 'add',items2test=(c(1:J)[-anchors])[which(c(1:J)[-anchors]%in%c(2,12:20))])#,anchor3,anchor4)])
  dif3.t=rbind(dif31[which(dif31$adj_pvals<0.05),],dif32[which(dif32$adj_pvals<0.05),])
  #power3[rep]=sum(c("V4" , "V5" , "V6", "V7","V8","V9","V12", "V13","V14","V15","V16","V17")%in%rownames(dif3.t))
  #tpI3[rep]=sum(c("V1","V2","V3","V10","V11","V18","V19","V20")%in%rownames(dif3.t))
  power3[rep]=sum(c("V4" , "V5" ,"V12", "V13")%in%rownames(dif3.t))
  tpI3[rep]=sum(c("V1","V2","V3", "V6", "V7","V8","V9","V10","V11","V14","V15","V16","V17","V18","V19","V20")%in%rownames(dif3.t))
  if ((power3[rep])==0){
    md.refit02 <-multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp02)[1:J]))
  } else {
    md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp02)[which(colnames(resp02)%in%rownames(dif3.t)==0)]))
  }
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  z23=(coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4:9,12:17)]
  difrec.mirt.fn3[rep,3]= mean(abs(z23[z23!=0]-0.6))
  z53=(coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[c(4:9)]
  z63=(coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[c(12:17)]
  difrec.mirt.fn3.a[rep,2]= mean(c(abs(z53[z53!=0]-0.6),abs(z63[z63!=0]-0.6)))
}

sum(power1)/(50*12)
sum(tpI1)/(50*8)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(na.omit(difrec.mirt.fn1))
colMeans(difrec.mirt.fn1.r)

sum(power2)/(50*12)
sum(tpI2)/(50*8)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

sum(power3)/(50*12)
sum(tpI3)/(50*8)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(difrec.mirt.fn3)

sd(power1/12)/7
sd(tpI1/8)/7
sd(power2/12)/7
sd(tpI2/8)/7
sd(power3/12)/7
sd(tpI3/8)/7


sum(power1)/(50*4)
sum(tpI1)/(50*16)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(na.omit(difrec.mirt.fn1))
colMeans(difrec.mirt.fn1.r)

sum(power2)/(50*4)
sum(tpI2)/(50*16)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

sum(power3)/(50*4)
sum(tpI3)/(50*16)
colMeans(bias.mirt3)
sqrt(colMeans((rmse.mirt3)^2))
colMeans(difrec.mirt.fn3)

sd(power1/4)/7
sd(tpI1/16)/7
sd(power2/4)/7
sd(tpI2/16)/7
sd(power3/4)/7
sd(tpI3/16)/7


write.csv(power1, file = "NALRTUnifpower1smallsample.csv")
write.csv(tpI1, file = "NALRTUniftpI1smallsample.csv")
write.csv(bias.mirt1, file = "NALRTUnifbias1smallsample.csv")
write.csv(rmse.mirt1, file = "NALRTUnifrmse1smallsample.csv")
write.csv(difrec.mirt.fn1, file = "NALRTUnifrec1smallsample.csv")
#write.csv(difrec.mirt.fn1.a, file = "NALRTNonUnifArec6lowcor.csv")
write.csv(power2, file = "NALRTUnifpower12smallsample.csv")
write.csv(tpI2, file = "NALRTUniftpI12smallsample.csv")
write.csv(difrec.mirt.fn2, file = "NALRTUnifrec12smallsamplecsv")
#write.csv(difrec.mirt.fn2.a, file = "NALRTNonUnifArec12smallsample.csv")
write.csv(power3, file = "NALRTUnifpower13smallsample.csv")
write.csv(tpI3, file = "NALRTUniftpI13smallsample.csv")
write.csv(difrec.mirt.fn3, file = "NALRTUnifrec13smallsample.csv")
#write.csv(difrec.mirt.fn3.a, file = "NALRTNonUnifArec13smallsample.csv")
write.csv(anchors.all, file = "NALRTUnifanchor1lowcor.csv")










