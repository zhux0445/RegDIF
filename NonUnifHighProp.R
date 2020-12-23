library(MASS)
library(mirt) 
mirtCluster(4)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(Rcpp)
library(RcppParallel)
sourceCpp("/Users/ruoyizhu/Documents/GitHub/mirt/matrix.cpp")
sourceCpp("/Users/zhux0445/Documents/GitHub/RegDIF/matrix.cpp")
setwd('/Users/zhux0445/Documents/GitHub/RegDIF_SimData')
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para4.csv",row.names = 1)
responses=read.csv("RESP6lowcor.csv",row.names = 1)

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

######################################################
###                                                ###
###                                                ###
###                50 replications                 ###
###                                                ###
###                                                ###
######################################################

reps=50
eta.5=numeric(reps)
Gammas.5=array(double(2*J*m*50),dim = c(2,2,J,50))
ADmat.5=array(double(J*3*reps),dim = c(J,3,reps)) #a has 2 columns, d has 1 column
biass.5=matrix(0,reps,3)
RMSEs.5=matrix(0,reps,3)

StartVals=read.csv("StartingValues6.csv",row.names = 1)
gra00=as.matrix(StartVals[,1:2])
rownames(gra00) <- c()
grd00=matrix(StartVals[,3],20,1)
grgamma00=array(0,dim=c(r,r,J))
grgamma00[c(1,2),1,3:11]=t(as.matrix(StartVals[3:11,4:5]))
grgamma00[c(1,2),2,12:20]=t(as.matrix(StartVals[12:20,4:5]))

# 2 dif per dim
mu100=c(0,0)
mu200=c(0,0.04)
mu300=c(0.01,0.01)
Sig100=matrix(c(1,0.8512375,0.8512375,1),2,2)
Sig200=matrix(c(1.19,1.05,1.05,1.3),2,2)
Sig300=matrix(c(0.91,0.87,0.87,1.15),2,2)

StartVals=read.csv("StartingValues6lowcor.csv")
gra00=as.matrix(StartVals[,1:2])
rownames(gra00) <- c()
grd00=matrix(StartVals[,3],20,1)
grgamma00=array(0,dim=c(r,r,J))
grgamma00[c(1,2),1,3:11]=t(as.matrix(StartVals[3:11,4:5]))
grgamma00[c(1,2),2,12:20]=t(as.matrix(StartVals[12:20,4:5]))

# 2 dif per dim
sigmat=read.csv("StartingValuesSig6lowcor.csv")
mu100=c(0,0)
mu200=c(0,0)
mu300=c(0,0)
Sig100=matrix(c(1,0.28,0.28,1),2,2)
Sig200=matrix(c(1.19,0.26,0.26,1.3),2,2)
Sig300=matrix(c(0.91,0.23,0.23,1.15),2,2)

#sim6 EM
for (rep in 2:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  r=2
  m=2
  eta.vec=seq(1,21,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  #Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=matrix(0,J,2),grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=T)
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Gammas[,,,k]=sim$Gamma
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  eta.5[rep]=eta.vec[kk]
  Gammas.5[,,,rep]=Gammas[,,,kk]
  ADmat.5[,,rep]=ADmat[,,kk]
  #biass.5[rep,]=biass[kk,]
  #RMSEs.5[rep,]=RMSEs[kk,]
  print(ADmat.5[,,rep])
  print(eta.5[rep])
  print(Gammas.5[,,,rep])
  #print(biass.5[rep,])
  #print(RMSEs.5[rep,])
  write.csv(eta.5[rep],file = paste("eta6Dec_",rep))
  write.csv(ADmat.5[,,rep],file = paste("ADmat6Dec_",rep))
  write.csv(rbind(t(rbind(Gammas.5[c(1,2),1,3:11,rep])),t(rbind(Gammas.5[c(1,2),2,12:20,rep]))),file = paste("Gamma6Dec_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta6Dec_",rep))
}

#sim6 lowcor EM
for (rep in 5:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  r=2
  m=2
  eta.vec=seq(3,21,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  #Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
 
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=matrix(0,J,2),grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=T)
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Gammas[,,,k]=sim$Gamma
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  eta.5[rep]=eta.vec[kk]
  Gammas.5[,,,rep]=Gammas[,,,kk]
  ADmat.5[,,rep]=ADmat[,,kk]
  #biass.5[rep,]=biass[kk,]
  #RMSEs.5[rep,]=RMSEs[kk,]
  print(ADmat.5[,,rep])
  print(eta.5[rep])
  print(Gammas.5[,,,rep])
  #print(biass.5[rep,])
  #print(RMSEs.5[rep,])
  write.csv(eta.5[rep],file = paste("eta6LowCor_",rep))
  write.csv(ADmat.5[,,rep],file = paste("ADmat6LowCor_",rep))
  write.csv(rbind(t(rbind(Gammas.5[c(1,2),1,3:11,rep])),t(rbind(Gammas.5[c(1,2),2,12:20,rep]))),file = paste("Gamma6LowCor_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta6LowCor_",rep))
}
