#MBR non-uniform DIF simulations
library(mirt) 
library(MASS)
library(mvtnorm)
mirtCluster(4)
library(ggplot2)

# Nonuniform DIF, 2D2PL 
#set.seed(1)
#A=runif(20,1.5,2.5)
A11=read.csv("Para1.csv",row.names = 1)[,1]
A12=A11;A13=A11
A12[4:9]=A12[4:9]-0.4;A13[4:9]=A13[4:9]-0.6
#A12[4:5]=A12[4:5]-0.4;A13[4:5]=A13[4:5]-0.6
A21=read.csv("Para1.csv",row.names = 1)[,2]
A22=A21;A23=A21
A22[12:17]=A22[12:17]-0.4;A23[12:17]=A23[12:17]-0.6
#A22[12:13]=A22[12:13]-0.4;A23[12:13]=A23[12:13]-0.6
D11=read.csv("Para1.csv",row.names = 1)[,7]
D12=D11;D13=D11
D12[c(4:9,12:17)]=D12[c(4:9,12:17)]+0.25
D13[c(4:9,12:17)]=D13[c(4:9,12:17)]+0.6
#D12[c(4:5,12:13)]=D12[c(4:5,12:13)]+0.25
#D13[c(4:5,12:13)]=D13[c(4:5,12:13)]+0.6

Amat1=cbind(A11,A21);Amat2=cbind(A12,A22);Amat3=cbind(A13,A23)
Dmat1=cbind(D11)
Dmat2=cbind(D12);Dmat3=cbind(D13)

#write.csv(cbind(Amat1,Amat2,Amat3,Dmat1,Dmat2,Dmat3), file = "Para4new.csv")
N1=N2=N3=1000
N=N1+N2+N3
set.seed(1)
#Theta=mvrnorm(N*50,c(0,0),matrix(c(1,0.25,0.25,1),2,2))
Theta=mvrnorm(N*50,c(0,0),matrix(c(1,0.85,0.85,1),2,2))
resp=matrix(0,N*50,J)
for (rep in 1:50){
  Theta1=Theta[((rep-1)*N+1):((rep-1)*N+N1),]
  Theta2=Theta[((rep-1)*N+N1+1):((rep-1)*N+N1+N2),]
  Theta3=Theta[((rep-1)*N+N1+N2+1):(rep*N),]
  datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
  datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
  datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
  resp[((rep-1)*N+1):(rep*N),]=rbind(datasetR,datasetF1,datasetF2)
}

write.csv(resp, file = "resp8new.csv")

# calculating wABC to check generated DIF magnitude
X1=seq(-4,4,by=0.25)
r=2
G=length(X1)^r
gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
X <- matrix(gh[idx,1],nrow(idx),r)
WXr=WXf1=dmvnorm(X,c(0,0),matrix(c(1,0.85,0.85,1),2,2))
G=length(X1)^r
m=2
J=20
gra=Amat1
grd=Dmat1
ygam1=Amat2-Amat1
#ygam1=Amat3-Amat1
ygamat1=ygam1%*%t(X)
ybeta1=Dmat2-Dmat1
#ybeta1=Dmat3-Dmat1
axmat=gra%*%t(X)
pstar1=pstar2=array(double(J*(m-1)*G),dim = c(J,(m-1),G))
p1=p2=array(double(J*m*G),dim = c(J,m,G))
for (g in 1:G)
{
  pstar1[,,g] = 1/(1+exp(-(grd+matrix(axmat[,g],J,1))))
  p1[,,g] = t(-diff(rbind(rep(1,J),t(pstar1[,,g]),rep(0,J)))) #(2)
}
for (g in 1:G)
{
  pstar2[,,g] = 1/(1+exp(-(grd+matrix(axmat[,g],J,1)+matrix(ygamat1[,g],J,1)+ybeta1))) #(3)
  p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
}

wABCr=wABCf1=numeric(J)
for (j in 1:J){
  wABCr[j]=sum(abs((p1[,2,]-p2[,2,])[j,])*WXr)
  wABCf1[j]=sum(abs((p1[,2,]-p2[,2,])[j,])*WXf1)
}
wABC1=  (wABCr+  wABCf1)/2 #for no impact condition, item 6,11,12,13,14 all have wABC between 0.3-0.9


library(MASS)
library(mirt) 
mirtCluster(4)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)
setwd('/Users/hyzhu27/Documents/GitHub/RegDIF_SimData')
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
setwd('/Users/zhux0445/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para3new.csv",row.names = 1)
responses=read.csv("RESP7new.csv",row.names = 1)

soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

J=20

N1=N2=N3=1000 
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


######################################################
###                                                ###
###                                                ###
###                50 replications                 ###
###                                                ###
###                                                ###
######################################################

reps=50
eta.2=numeric(reps)
Gammas.2=array(double(2*J*m*50),dim = c(2,2,J,50))
Betas.2=array(double(J*2*reps),dim = c(J,2,reps))
ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) #a has 2 columns, d has 1 column
biass.2=matrix(0,reps,3)
RMSEs.2=matrix(0,reps,3)

s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'

starting7new=matrix(0,20,7)
# starting value method 1

for(i in 1:20){
  md.cons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[-i]))
  starting5new[i,1:3]=coef(md.cons0,simplify=T)$G1$items[i,1:3]
  starting5new[i,4]=sum(coef(md.cons0,simplify=T)$G2$items[i,1:2]-coef(md.cons0,simplify=T)$G1$items[i,1:2])
  starting5new[i,5]=sum(coef(md.cons0,simplify=T)$G3$items[i,1:2]-coef(md.cons0,simplify=T)$G1$items[i,1:2])
  starting5new[i,6]=coef(md.cons0,simplify=T)$G2$items[i,3]-coef(md.cons0,simplify=T)$G1$items[i,3]
  starting5new[i,7]=coef(md.cons0,simplify=T)$G3$items[i,3]-coef(md.cons0,simplify=T)$G1$items[i,3]
}

# starting value method 2
md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE)
coef(md.noncons0,simplify=T)
mean((coef(md.noncons0,simplify=T)$G1$items[,1]-Amat1[,1])^2)
mean((gra00[,1]-Amat1[,1])^2)

mean((coef(md.noncons0,simplify=T)$G1$items[,3]-Dmat1[,1])^2)
mean((grd00[,1]-Dmat1[,1])^2)

cbind(coef(md.noncons0,simplify=T)$G2$items[,3]-coef(md.noncons0,simplify=T)$G1$items[,3],coef(md.noncons0,simplify=T)$G3$items[,3]-coef(md.noncons0,simplify=T)$G1$items[,3])

cbind(coef(md.noncons0,simplify=T)$G2$items[,1]-coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G2$items[,2]-coef(md.noncons0,simplify=T)$G1$items[,2])
cbind(coef(md.noncons0,simplify=T)$G3$items[,1]-coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G3$items[,2]-coef(md.noncons0,simplify=T)$G1$items[,2])
starting5new=cbind(coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G1$items[,2],coef(md.noncons0,simplify=T)$G1$items[,3],
                   rowSums(cbind(coef(md.noncons0,simplify=T)$G2$items[,1]-coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G2$items[,2]-coef(md.noncons0,simplify=T)$G1$items[,2])),
rowSums(cbind(coef(md.noncons0,simplify=T)$G3$items[,1]-coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G3$items[,2]-coef(md.noncons0,simplify=T)$G1$items[,2])),
cbind(coef(md.noncons0,simplify=T)$G2$items[,3]-coef(md.noncons0,simplify=T)$G1$items[,3],coef(md.noncons0,simplify=T)$G3$items[,3]-coef(md.noncons0,simplify=T)$G1$items[,3]))

#write.csv(cbind(gra00,grd00,grbeta00),file = "StartingValues2.csv")
gra00=as.matrix(starting5new[,1:2])
rownames(gra00) <- c()
grd00=matrix(starting5new[,3],20,1)
grgamma00=array(0,dim=c((y-1),r,J))
grgamma00[1,1,1]=starting5new[1,4];grgamma00[1,2,2]=starting5new[2,4];grgamma00[1,1,3:11]=starting5new[3:11,4];grgamma00[1,2,12:20]=starting5new[12:20,4]
grgamma00[2,1,1]=starting5new[1,5];grgamma00[2,2,2]=starting5new[2,5];grgamma00[2,1,3:11]=starting5new[3:11,5];grgamma00[2,2,12:20]=starting5new[12:20,5]
grbeta00=as.matrix(starting5new[,6:7])

# 6 dif per dim
mu100=c(0,0)
mu200=c(0.1472920, -0.0157996)
mu300=c(0.1100240, -0.1232052)
Sig100=matrix(c(1,0.2753316,0.2753316,1),2,2)
Sig200=matrix(c(1.3259608,0.3145355,0.3145355,1.1363796),2,2)
Sig300=matrix(c(1.2270710,0.2503095,0.2503095,1.0718629),2,2)

mu100=c(0,0)
mu200=c(0,0)
mu300=c(0,0)
Sig100=Sig200=Sig300=matrix(c(1,0.8512375,0.8512375,1),2,2)
Sig200=matrix(c(1,0.846,0.846,1),2,2)
Sig300=matrix(c(1,0.858,0.858,1),2,2)



#sim1 EM
for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  r=2
  m=2
  y=3
  eta.vec=seq(21,45,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    print(proc.time() - ptm)
    bics[k]=sim$bic
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    Gammas[,,,k]=sim$Gamma
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  Gammas.2[,,,rep]=Gammas[,,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(Gammas.2[,,,rep])
  write.csv(eta.2[rep],file = paste("NAeta7_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat7_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta7_",rep))
  write.csv(rbind(t(rbind(Gammas.2[c(1,2),1,c(1,3:11),rep])),t(rbind(Gammas.2[c(1,2),2,c(2,12:20),rep]))),file = paste("NAGamma7_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta7_",rep))
}

#sim1 EMM
for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  r=2
  m=2
  y=3
  eta.vec=seq(21,45,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_EMM_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    print(proc.time() - ptm)
    bics[k]=sim$bic
    Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  Gammas.2[,,,rep]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("NAeta7EMM_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat7EMM_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta7EMM_",rep))
  write.csv(rbind(t(rbind(Gammas.2[c(1,2),1,c(1,3:11),rep])),t(rbind(Gammas.2[c(1,2),2,c(2,12:20),rep]))),file = paste("NAGamma7EMM_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta7EMM_",rep))
}

#sim1 adaptive
for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  r=2
  m=2
  y=3
  eta.vec=seq(1,25,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_Adaptive_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,lam=1,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    
    print(proc.time() - ptm)
    bics[k]=sim$bic
    Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  Gammas.2[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  write.csv(eta.2[rep],file = paste("NAeta5Adapt_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat5Adapt_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta5Adapt_",rep))
  write.csv(rbind(t(rbind(Gammas.2[c(1,2),1,c(1,3:11),rep])),t(rbind(Gammas.2[c(1,2),2,c(2,12:20),rep]))),file = paste("NAGamma5Adapt_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta5Adapt_",rep))
}
