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
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("RESP3.csv",row.names = 1)

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
#Gammas.1=array(double(2*J*m*50),dim = c(2,2,J,50))
Betas.2=array(double(J*2*reps),dim = c(J,2,reps))
ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) #a has 2 columns, d has 1 column
biass.2=matrix(0,reps,3)
RMSEs.2=matrix(0,reps,3)

#write.csv(cbind(gra00,grd00,grbeta00),file = "StartingValues2.csv")
StartVals=read.csv("StartingValues1.csv",row.names = 1)
gra00=as.matrix(StartVals[,1:2])
rownames(gra00) <- c()
grd00=matrix(StartVals[,3],20,1)
grbeta00=as.matrix(StartVals[,4:5])
grbeta00[1,]=c(-0.0481844253,0.02196438);grbeta00[2,]=c(0.1594530242,0.15120530)
rownames(grbeta00) <- c()
colnames(grbeta00) <- c()

# 6 dif per dim
mu100=c(0,0)
mu200=c(0.1472920, -0.0157996)
mu300=c(0.1100240, -0.1232052)
Sig100=matrix(c(1,0.2753316,0.2753316,1),2,2)
Sig200=matrix(c(1.3259608,0.3145355,0.3145355,1.1363796),2,2)
Sig300=matrix(c(1.2270710,0.2503095,0.2503095,1.0718629),2,2)

StartVals=read.csv("StartingValues2.csv",row.names = 1)
gra00=as.matrix(StartVals[,1:2])
rownames(gra00) <- c()
grd00=matrix(StartVals[,3],20,1)
grbeta00=as.matrix(StartVals[,4:5])
grbeta00[1,]=c(-0.0481844253,0.02196438);grbeta00[2,]=c(0.1594530242,0.15120530)
rownames(grbeta00) <- c()
colnames(grbeta00) <- c()

StartVals=read.csv("StartingValues2.csv",row.names = 1)
gra00=as.matrix(StartVals[,1:2])
rownames(gra00) <- c()
grd00=matrix(StartVals[,3],20,1)
grbeta00=as.matrix(StartVals[,4:5])
rownames(grbeta00) <- c()
colnames(grbeta00) <- c()
mu100=c(0,0)
mu200=c(-0.1022104,0.1383983)
mu300=c(0.03270989,0.06991659)
Sig100=matrix(c(1,0.8512375,0.8512375,1),2,2)
Sig200=matrix(c(0.9879547,0.8953391,0.8953391,1.0742965),2,2)
Sig300=matrix(c(0.8755486,0.8193335,0.8193335,1.0120597),2,2)


#########################
# 6 dif per dim
mu100=c(0,0)
mu200=c(0,0)
mu300=c(0,0)
Sig100=matrix(c(1,0.2892655,0.2892655,1),2,2)
Sig200=matrix(c(1.0518386,0.2419183,0.2419183,1.0355795),2,2)
Sig300=matrix(c(0.9386327,0.2428997,0.2428997,0.9842235),2,2)

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
  eta.vec=seq(11,35,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("NAeta3_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat3_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta3_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta3_",rep))
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
  eta.vec=seq(11,35,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_EMM_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("NAeta1EMM_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat1EMM_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta1EMM_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta1EMM_",rep))
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
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_Adaptive_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,lam=1,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
      
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  write.csv(eta.2[rep],file = paste("NAeta1Adapt_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat1Adapt_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta1Adapt_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta1Adapt_",rep))
}

#sim2 EM
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
  eta.vec=seq(11,35,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("NAeta2_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat2_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta2_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta2_",rep))
}

#sim2 EMM
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
  eta.vec=seq(11,35,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_EMM_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("NAeta1EMM_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat1EMMLowCor_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta1EMMLowCor_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta1EMMLowCor_",rep))
}

#sim2 adaptive
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
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_Adaptive_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,lam=1,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
    
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  write.csv(eta.2[rep],file = paste("NAeta1Adapt_",rep))
  write.csv(ADmat.2[,,rep],file = paste("NAADmat1Adapt_",rep))
  write.csv(Betas.2[,,rep],file = paste("NABeta1Adapt_",rep))
  write.csv(theta.dist[,,kk],file = paste("NAtheta1Adapt_",rep))
}

# sim3 Lowcor
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
  y=3
  eta.vec=seq(15,75,5)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
    print(proc.time() - ptm) 
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  write.csv(eta.2[rep],file = paste("eta3LowCor_",rep))
  write.csv(ADmat.2[,,rep],file = paste("ADmat3LowCor_",rep))
  write.csv(Betas.2[,,rep],file = paste("Beta3LowCor_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta3LowCor_",rep))
}

# sim3 Lowcor EMM
for (rep in 31:reps){
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
  eta.vec=seq(15,75,5)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_EMM_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
    print(proc.time() - ptm) 
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  write.csv(eta.2[rep],file = paste("eta3EMM_",rep))
  write.csv(ADmat.2[,,rep],file = paste("ADmat3EMM_",rep))
  write.csv(Betas.2[,,rep],file = paste("Beta3EMM_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta3EMM_",rep))
}

#sim3 lower adaptive
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
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_Adaptive_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,lam=1,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
    
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  write.csv(eta.2[rep],file = paste("eta3AdaptLowCor_",rep))
  write.csv(ADmat.2[,,rep],file = paste("ADmat3AdaptLowCor_",rep))
  write.csv(Betas.2[,,rep],file = paste("Beta3AdaptLowCor_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta3AdaptLowCor_",rep))
  
}




#sim4 EM
for (rep in 41:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  r=2
  m=2
  y=3
  eta.vec=seq(1,31,3)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(500,500,500),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("eta2EMHighCor_",rep))
  write.csv(ADmat.2[,,rep],file = paste("ADmat2EMHighCor_",rep))
  write.csv(Betas.2[,,rep],file = paste("Beta2EMHighCor_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta2EMHighCor_",rep))
}


#sim4 lower EM
for (rep in 24:reps){
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
  eta.vec=seq(10,37,3)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  theta.dist=array(double(2*9*length(eta.vec)),dim=c(9,2,length(eta.vec)))
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(1000,1000,1000),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=array(0,dim=c((y-1),r,J)),Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300),NonUniform=F)
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    theta.dist[,,k]=rbind(sim$mean1,sim$mean2,sim$mean3,sim$Corr1,sim$Corr2,sim$Corr3)
  }
  
  kk=which.min(bics)
  
  eta.2[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.2[,,rep]=ADmat[,,kk]
  Betas.2[,,rep]=Betas[,,kk]
  biass.2[rep,]=biass[kk,]
  RMSEs.2[rep,]=RMSEs[kk,]
  print(ADmat.2[,,rep])
  print(eta.2[rep])
  print(Betas.2[,,rep])
  print(biass.2[rep,])
  print(RMSEs.2[rep,])
  write.csv(eta.2[rep],file = paste("eta4EM_",rep))
  write.csv(ADmat.2[,,rep],file = paste("ADmat4EM_",rep))
  write.csv(Betas.2[,,rep],file = paste("Beta4EM_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta4EM_",rep))
}















