# DIF identification
# Regularized MGRM

library(MASS)
library(mirt)
# num of items is fixed to 20
J=20

# Dataset #1 (2 Uniform DIF items per scale)
# item (4,5,12,13) have DIF and item (1,11) are anchor (DIF free)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A13=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A23=A21
D11=rnorm(J,0,1)
D12=D13=D11
# focal group has more extreme intercepts
D12[4]=D11[4]+ 0.5
D12[5]=D11[5]+ 0.5
D13[4]=D11[4]+ 1
D13[5]=D11[5]+ 1
D12[12]=D11[12]+ 0.5
D12[13]=D11[13]+ 0.5
D13[12]=D11[12]+ 1
D13[13]=D11[13]+ 1
Amat01=cbind(A11,A21)
Amat1=Amat01
Amat1[2,]=Amat01[11,]
Amat1[11,]=Amat01[2,]
Amat02=cbind(A12,A22)
Amat2=Amat02
Amat2[2,]=Amat02[11,]
Amat2[11,]=Amat02[2,]
Amat03=cbind(A13,A23)
Amat3=Amat03
Amat3[2,]=Amat03[11,]
Amat3[11,]=Amat03[2,]
Dmat1=cbind(D11)
Dmat2=cbind(D12)
Dmat3=cbind(D13)

# Dataset #2 (4 Uniform DIF items per scale)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A21
D11=rnorm(J,0,1)
D12=D11
# focal group has more extreme intercepts
D12[3]=D11[3]+ 0.25
D12[4]=D11[4]+ 0.25 
D12[5]=D11[5]+ 0.5
D12[6]=D11[6]+ 0.5
D12[11]=D11[11]+ 0.25
D12[12]=D11[12]+ 0.25
D12[13]=D11[13]+ 0.5
D12[14]=D11[14]+ 0.5
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Dmat1=cbind(D11)
Dmat2=cbind(D12)

# Dataset #3 (2 Non-uniform DIF items per scale)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A21
D11=rnorm(J,0,1)
D12=D11
# focal group has smaller slopes and more extreme intercepts
A12[4]=A11[4]- 0.3
A12[5]=A11[5]- 0.6
A22[12]=A21[12]- 0.3
A22[13]=A21[13]- 0.6
#D12[4]=D11[4]+ 0.25
#D12[5]=D11[5]+ 0.5
#D12[12]=D11[12]+ 0.25
#D12[13]=D11[13]+ 0.5
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Dmat1=cbind(D11)
Dmat2=cbind(D12)

# Dataset #4 (4 Non-uniform DIF items per scale)
set.seed(1)
A11=c(runif(J/2,1.5,2.5),rep(0,J/2))
A12=A11
A21=c(rep(0,J/2),runif(J/2,1.5,2.5))
A22=A21
D11=rnorm(J,0,1)
D12=D11
# focal group has smaller slopes and more extreme intercepts
A12[3]=A11[3]-0.3
A12[4]=A11[4]-0.3
A12[5]=A11[5]-0.6
A12[6]=A11[6]-0.6
A22[11]=A21[11]- 0.3
A22[12]=A21[12]- 0.3
A22[13]=A21[13]- 0.6
A22[14]=A21[14]- 0.6
D12[3]=D11[3]+ 0.25
D12[4]=D11[4]+ 0.25
D12[5]=D11[5]+ 0.5
D12[6]=D11[6]+ 0.5
D12[11]=D11[11]+ 0.25
D12[12]=D11[12]+ 0.25
D12[13]=D11[13]+ 0.5
D12[14]=D11[14]+ 0.5
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Dmat1=cbind(D11)
Dmat2=cbind(D12)


N1=N2=N3=500
N=N1+N2+N3
Group=c(rep('G1', N1), rep('G2', N2),rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
y1=c(0,0)
y2=matrix(c(1,0),1,2)
y3=matrix(c(0,1),1,2)
y=matrix(0,N1+N2+N3,2)
y[1:N1,]=y1
for (i in (N1+1):(N1+N2)){
  y[i,]=y2
}
for (i in (N1+N2+1):(N1+N2+N3)){
  y[i,]=y3
}

set.seed(200)
Theta=mvrnorm(n=N1+N2+N3,mu=c(0,0),Sigma = matrix(c(1,0.85,0.85,1),2,2))
Theta1=Theta[1:N1,]
Theta2=Theta[(N1+1):(N1+N2),]
Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
resp=rbind(datasetR,datasetF1,datasetF2)

# CTT starting values
total.sco = rowSums(resp)
sy=sd(total.sco)
totol.mat=matrix(0,2,J)
num.mat=matrix(0,2,J)

for (j in 1:J){
  totol.mat[1,j] = mean(total.sco[which(resp[,j]==0)])
  totol.mat[2,j] = mean(total.sco[which(resp[,j]==1)])
  num.mat[1,j] = length(total.sco[which(resp[,j]==0)])
  num.mat[2,j] = length(total.sco[which(resp[,j]==1)])
}

rpb=(num.mat[1,]*num.mat[2,]/(dnorm(num.mat[2,]/N)*N^2))*(totol.mat[2,]-totol.mat[1,])/sy
a0=rpb/sqrt(1-rpb^2) 
azero=matrix(0,J,r)
azero[1,1]=a0[1]
azero[2,2]=a0[2]
azero[3:11,1]=a0[3:11]
azero[12:20,2]=a0[12:20]
b0=-(1/qnorm(num.mat[2,]/N))/rpb 
d0=-a0*b0 #d looks incorrect

r=2
m=2
eta.vec=seq(30,40,2)
bics=rep(0,length(eta.vec))
Betas=array(double(J*m*length(eta.vec)),dim = c(J,m,length(eta.vec)))
biass=matrix(0,length(eta.vec),3)
RMSEs=matrix(0,length(eta.vec),3)
for (k in 1:length(eta.vec))
{
  eta=eta.vec[k]
  sim=ipest1(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=gra00,grd00=grd00,grbeta00=grbeta00,Sig00=Sig00)
  bics[k]=sim$bic
  Betas[,,k]=sim$Beta
  biass[k,]=sim$bias
  RMSEs[k,]=sim$RMSE
}

kk=which.min(bics)
eta.vec[kk]
Betas[,,kk]
biass[kk,]
RMSEs[kk,]

colSums(gra-Amat1)/10
sqrt(colSums((gra-Amat1)^2)/10)
colMeans(grd-Dmat1)
sqrt(colMeans((grd-Dmat1)^2))
Sig

colMeans(grd+grbeta[,1]-Dmat2)
sqrt(colMeans((grd+grbeta[,1]-Dmat2)^2))

colMeans(grd+grbeta[,2]-Dmat3)
sqrt(colMeans((grd+grbeta[,2]-Dmat3)^2))







#########################################
#####                               #####
#####  Comparing with mirt LRT add  #####
#####                               #####
#########################################

mirt.p.mat1=matrix(0,16,18) # 16 is the number of replications
bias.mirt1=matrix(0,16,3)
rmse.mirt1=matrix(0,16,3)
difrec.mirt.fn1=matrix(0,16,3) # DIF magnitude recovery with false negative included
mirt.p.mat2=matrix(0,16,18) # 16 is the number of replications
bias.mirt2=matrix(0,16,3)
rmse.mirt2=matrix(0,16,3)
difrec.mirt.fn2=matrix(0,16,3) # 
mirt.p.mat3=matrix(0,16,18) # 16 is the number of replications
bias.mirt3=matrix(0,16,3)
rmse.mirt3=matrix(0,16,3)
difrec.mirt.fn3=matrix(0,16,3) 

for (i in 1:16){
  seed=i*100
  set.seed(seed)
  Theta=mvrnorm(n=N1+N2+N3,mu=c(0,0),Sigma = matrix(c(1,0.85,0.85,1),2,2))
  Theta1=Theta[1:N1,]
  Theta2=Theta[(N1+1):(N1+N2),]
  Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
  datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
  datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
  datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
  resp=rbind(datasetR,datasetF1,datasetF2)
  resp01=rbind(datasetR,datasetF1)
  resp02=rbind(datasetR,datasetF2)
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:r]))
  mirt.p.mat1[(i),]=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-(which(mirt.p.mat1[i,]<0.05)+2)]))
  bias.mirt1[i,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[i,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[i,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[i,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[i,2:3]= c(mean(abs((coef(md.refit0,simplify=T)$G2$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit0,simplify=T)$G3$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  difrec.mirt.fn1[i,1]=mean(c(difrec.mirt.fn[i,2],difrec.mirt.fn[i,3]))
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:r]))
  mirt.p.mat2[(i),]=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[-(which(mirt.p.mat2[i,]<0.05)+2)]))
  bias.mirt2[i,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[i,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[i,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[i,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[i,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:r]))
  mirt.p.mat3[(i),]=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[-(which(mirt.p.mat3[i,]<0.05)+2)]))
  bias.mirt3[i,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[i,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[i,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[i,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[i,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1))
  }
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(16*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(16*14)
colMeans(bias.mirt1)
colMeans(rmse.mirt1)
colMeans(difrec.mirt.fn1)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(16*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(16*14)
colMeans(bias.mirt2)
colMeans(rmse.mirt2)
colMeans(difrec.mirt.fn2)

sum(mirt.p.mat3[,c(2,3,10,11)]<0.05)/(16*4)
sum(mirt.p.mat3[,-c(2,3,10,11)]<0.05)/(16*14)
colMeans(bias.mirt3)
colMeans(rmse.mirt3)
colMeans(difrec.mirt.fn3)

#########################################
#####                               #####
#####  Comparing with mirt Wald     #####
#####                               #####
#########################################

mirt.p.mat1w=matrix(0,16,18) # 16 is the number of replications
bias.mirt1w=matrix(0,16,3)
rmse.mirt1w=matrix(0,16,3)
difrec.mirt.fn1w=matrix(0,16,3) # DIF magnitude recovery with false negative included
mirt.p.mat2w=matrix(0,16,18) # 16 is the number of replications
bias.mirt2w=matrix(0,16,3)
rmse.mirt2w=matrix(0,16,3)
difrec.mirt.fn2w=matrix(0,16,3) # 
mirt.p.mat3w=matrix(0,16,18) # 16 is the number of replications
bias.mirt3w=matrix(0,16,3)
rmse.mirt3w=matrix(0,16,3)
difrec.mirt.fn3w=matrix(0,16,3) 

for (i in 1:16){
  seed=i*100
  set.seed(seed)
  Theta=mvrnorm(n=N1+N2+N3,mu=c(0,0),Sigma = matrix(c(1,0.85,0.85,1),2,2))
  Theta1=Theta[1:N1,]
  Theta2=Theta[(N1+1):(N1+N2),]
  Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
  datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
  datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
  datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
  resp=rbind(datasetR,datasetF1,datasetF2)
  resp01=rbind(datasetR,datasetF1)
  resp02=rbind(datasetR,datasetF2)
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:r]))
  mirt.p.mat1w[(i),]=DIF(md.noncons0, which.par = c('d'), Wald = TRUE, p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-(which(mirt.p.mat1w[i,]<0.05)+2)]))
  bias.mirt1w[i,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1w[i,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1w[i,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1w[i,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1w[i,2:3]= c(mean(abs((coef(md.refit0,simplify=T)$G2$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit0,simplify=T)$G3$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  difrec.mirt.fn1w[i,1]=mean(c(difrec.mirt.fn[i,2],difrec.mirt.fn[i,3]))
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:r]))
  mirt.p.mat2w[(i),]=DIF(md.noncons01, which.par = c('d'), Wald = TRUE, p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[-(which(mirt.p.mat2w[i,]<0.05)+2)]))
  bias.mirt2w[i,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2w[i,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2w[i,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2w[i,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2w[i,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:r]))
  mirt.p.mat3w[(i),]=DIF(md.noncons02, which.par = c('d'), Wald = TRUE, p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[-(which(mirt.p.mat3w[i,]<0.05)+2)]))
  bias.mirt3w[i,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3w[i,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3w[i,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3w[i,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3w[i,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1))
}
sum(mirt.p.mat1w[,c(2,3,10,11)]<0.05)/(16*4)
sum(mirt.p.mat1w[,-c(2,3,10,11)]<0.05)/(16*14)
colMeans(bias.mirt1w)
colMeans(rmse.mirt1w)
colMeans(difrec.mirt.fn1w)

sum(mirt.p.mat2w[,c(2,3,10,11)]<0.05)/(16*4)
sum(mirt.p.mat2w[,-c(2,3,10,11)]<0.05)/(16*14)
colMeans(bias.mirt2w)
colMeans(rmse.mirt2w)
colMeans(difrec.mirt.fn2w)

sum(mirt.p.mat3w[,c(2,3,10,11)]<0.05)/(16*4)
sum(mirt.p.mat3w[,-c(2,3,10,11)]<0.05)/(16*14)
colMeans(bias.mirt3w)
colMeans(rmse.mirt3w)
colMeans(difrec.mirt.fn3w)

# cannot save p-values appropriatly.
#ref vs focal1
md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('slopes',colnames(resp01)[-(which(mirt.p.mat1[i,]<0.05)+2)]))
mirt.p.mat2[(i),]=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
#ref vs focal2
md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('slopes',colnames(resp02)[-(which(mirt.p.mat1[i,]<0.05)+2)]))
mirt.p.mat3[(i),]=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]


itemname3 <- colnames(dat)
refmod.noncon3 <- multipleGroup(dat, model3, group = Group, SE=TRUE,invariance=c('free_means', 'free_var', itemname3[c(1,21)]))
# loop over items 
estmodels <- vector('list', ncol(dat)-2) # minus 2 anchors
for(i in 1:(ncol(dat)-2))
  estmodels[[i]] <- multipleGroup(dat, model3, group = Group, verbose = FALSE,calcNull=FALSE,invariance=c('free_means', 'free_var', itemname3[c(1,i+1,21)]))
(anovas <- lapply(estmodels, anova, object2=refmod.noncon3, verbose=FALSE)) 
p <- do.call(rbind, lapply(anovas, function(x) x[2, 'p']))
p.adjust(p, method = 'BH')
anchoritem3=which(p.adjust(p, method = 'BH')>0.05)+1
refitmodel3=multipleGroup(dat, model3, group = Group, SE=TRUE,invariance=c('free_means', 'free_var', itemname3[c(1,anchoritem3,21)]))
coef(refitmodel3)
