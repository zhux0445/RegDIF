# Theta estimation
library(MASS)
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)

MGRM_est <- function(resp,params,ncat,p,constr,pr_mu, invcov){
  # Preparations
  L = length(resp) #get the current test length
  if(length(ncat)==1){ncat=rep(ncat,L)}
  if(length(ncat)!=L){
    stop("Incorrect response categories")
  } # check response categories
  #if(ncol(params)!=p+max(ncat)-1){
  #stop("Incorrect Item Parameters")
  #} # check item parameters
  a = params[,c(1:p)] #extract first p rows of the matrix params, 
  # the discriminating parameters, a is a matrix of size L by p
  b = params[,c((p+1):ncol(params))] # extract the difficulty parameter matrix
  b = cbind(b,rep(-1e8,L)) # add one column to the right of b (to vectorize the likelihood computation)
  b = cbind(rep(1e8,L),b) # add one column to the left of b (to vectorize the likelihood computation)
  ## MAP Estimation 
  ## this is the same to MLE except the prior information is not 0 (the first step in computing 
  ## the objective function)
  idx1 = cbind(c(1:L),resp+1)
  idx2 = idx1
  idx2[,2] = idx2[,2]+1
  b1 = matrix(rep(b[idx1],p),L,1)
  b2 = matrix(rep(b[idx2],p),L,1)
  obj <- function(x){
    #temp = 0 # MLE
    temp = -0.5*crossprod(x-pr_mu,invcov)%*%(x-pr_mu) #prior dist of theta for MAP
    for(j in 1:L){
      temp=temp+log(1/(1+exp(-(a[j,]%*%x+b1[j,])))-
                      1/(1+exp(-(a[j,]%*%x+b2[j,]))))
    }
    return(-temp)
  }
  temp = optim(rep(0,p),obj,
               lower=rep(constr[1],p),upper=rep(constr[2],p),
               method="L-BFGS-B",hessian=TRUE)
  theta_hat = temp$par
  fimat = solve(temp$hessian)
  std_error = sqrt(diag(fimat))
  detfi = det(temp$hessian)
  return(list(theta_hat=theta_hat,std_error = std_error, detfi=detfi))
}

## sim 1
J=20
N1=N2=N3=1000
N=N1+N2+N3
set.seed(100)
Theta=mvrnorm(n=N*50,mu=c(0,0),Sigma = matrix(c(1,0.25,0.25,1),2,2))
# Check whether theta's are same as those in RegDIF simulations.
Theta1=Theta[1:N1,]
Theta2=Theta[(N1+1):(N1+N2),]
Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para4.csv",row.names = 1)
Amat1=as.matrix(params[,1:2])
rownames(Amat1) <- c()
Amat2=as.matrix(params[,3:4])
rownames(Amat2) <- c()
Amat3=as.matrix(params[,5:6])
rownames(Amat3) <- c()
Dmat1=matrix(params[,7],J,(m-1))
Dmat2=matrix(params[,8],J,(m-1))
Dmat3=matrix(params[,9],J,(m-1))
reps=50
resp.all=matrix(0,N*reps,J)
for (rep in 1:reps){
  Theta1=Theta[((rep-1)*N+1):((rep-1)*N+N1),]
  Theta2=Theta[((rep-1)*N+N1+1):((rep-1)*N+N1+N2),]
  Theta3=Theta[((rep-1)*N+N1+N2+1):((rep-1)*N+N1+N2+N3),]
  datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
  datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
  datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
  resp.all[((rep-1)*N+1):((rep-1)*N+N),]=rbind(datasetR,datasetF1,datasetF2)
}
write.csv(resp.all,file = "RESP8lowcor.csv")
responses=as.matrix(read.csv("RESP3.csv",row.names = 1))
sum(resp-responses[1:3000,])
colnames(responses)=c()
rownames(responses)=c()

reps=50
theta.est.mat2=Theta
theta.se.mat2=Theta
para.mat.gp1=para.mat.gp2=para.mat.gp3=matrix(0,J,3)
for (rep in 1:50){
  resp1=responses[((rep-1)*N+1):((rep-1)*N+N1),]
  resp2=responses[((rep-1)*N+N1+1):((rep-1)*N+N1+N2),]
  resp3=responses[((rep-1)*N+N1+N2+1):((rep-1)*N+N1+N2+N3),]
  para.mat.gp1=para.mat.gp2=para.mat.gp3=as.matrix(read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/ADmat3_",rep),row.names = 1))
  colnames(para.mat.gp1)=c();colnames(para.mat.gp2)=c();colnames(para.mat.gp3)=c()
  rownames(para.mat.gp1)=c();rownames(para.mat.gp2)=c();rownames(para.mat.gp3)=c()
  para.mat.gp2[,3]=para.mat.gp1[,3]+read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[,1]
  para.mat.gp3[,3]=para.mat.gp1[,3]+read.csv(paste("/Users/ruoyizhu/Desktop/condition 3_all replications/Beta3_",rep),row.names = 1)[,2]
  for (i in 1:N1){
    #theta.esti=MGRM_est(resp=resp1[i,],params=para.mat.gp1,ncat=2,p=2,constr = c(-5,5)) # mle
    theta.esti=MGRM_est(resp=resp1[i,],params=para.mat.gp1,ncat=2,p=2,constr = c(-5,5),pr_mu=c(0,0), invcov=matrix(c(2.77,-2.222,-2.222,2.778),2,2)) # map
    theta.est.mat2[(rep-1)*N+i,]=theta.esti$theta_hat
    theta.se.mat2[(rep-1)*N+i,]=theta.esti$std_error
  }
  for (i in 1:N2){
    #theta.esti=MGRM_est(resp=resp2[i,],params=para.mat.gp2,ncat=2,p=2,constr = c(-5,5))# mle
    theta.esti=MGRM_est(resp=resp2[i,],params=para.mat.gp2,ncat=2,p=2,constr = c(-5,5),pr_mu=c(0,0), invcov=matrix(c(2.77,-2.222,-2.222,2.778),2,2))# map
    theta.est.mat2[(rep-1)*N+N1+i,]=theta.esti$theta_hat
    theta.se.mat2[(rep-1)*N+N1+i,]=theta.esti$std_error
  }
  for (i in 1:N3){
    #theta.esti=MGRM_est(resp=resp3[i,],params=para.mat.gp3,ncat=2,p=2,constr = c(-5,5)) #mle
    theta.esti=MGRM_est(resp=resp3[i,],params=para.mat.gp3,ncat=2,p=2,constr = c(-5,5),pr_mu=c(0,0), invcov=matrix(c(2.77,-2.222,-2.222,2.778),2,2)) #map
    theta.est.mat2[(rep-1)*N+N1+N2+i,]=theta.esti$theta_hat
    theta.se.mat2[(rep-1)*N+N1+N2+i,]=theta.esti$std_error
  }
}
colMeans(theta.est.mat2-Theta)
sqrt(colMeans((theta.est.mat2-Theta)^2))
colMeans(theta.se.mat2)

mirt.p.mat1=read.csv("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Sim7_LRTpvs1.csv",row.names = 1)
sum(mirt.p.mat1[1:8,c(2,3,10,11)]<0.05)/(8*4)
