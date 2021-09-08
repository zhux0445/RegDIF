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


params=read.csv("Para3new.csv",row.names = 1)
responses=read.csv("RESP5new.csv",row.names = 1)

J=20 #test length
N1=N2=N3=500 # sample size each group
N=N1+N2+N3

m=2 # number of response categories
r=2 # number of latent factors
y=3 # number of covariate groups
y1=c(0,0) # group indicator 
y2=c(1,0)
y3=c(0,1)

Amat1=params[,1:2] # set matrices for generated item parameters
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



reps=50 # total number of replications
rep=1 # take the first replication as example in analysis below 
resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
eta.2=numeric(reps)
Gammas.2=array(double(2*J*m*50),dim = c(2,2,J,50))
Betas.2=array(double(J*2*reps),dim = c(J,2,reps))
ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) #a has 2 columns, d has 1 column
biass.2=matrix(0,reps,3)
RMSEs.2=matrix(0,reps,3)

s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'

# calculate starting value using multiplegroup
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE)
starting5new=cbind(coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G1$items[,2],coef(md.noncons0,simplify=T)$G1$items[,3],
                   rowSums(cbind(coef(md.noncons0,simplify=T)$G2$items[,1]-coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G2$items[,2]-coef(md.noncons0,simplify=T)$G1$items[,2])),
                   rowSums(cbind(coef(md.noncons0,simplify=T)$G3$items[,1]-coef(md.noncons0,simplify=T)$G1$items[,1],coef(md.noncons0,simplify=T)$G3$items[,2]-coef(md.noncons0,simplify=T)$G1$items[,2])),
                   cbind(coef(md.noncons0,simplify=T)$G2$items[,3]-coef(md.noncons0,simplify=T)$G1$items[,3],coef(md.noncons0,simplify=T)$G3$items[,3]-coef(md.noncons0,simplify=T)$G1$items[,3]))

# starting values estimated by multiplegroup function
gra00=as.matrix(starting5new[,1:2])
rownames(gra00) <- c()
grd00=matrix(starting5new[,3],20,1)
grgamma00=array(0,dim=c((y-1),r,J))
grgamma00[1,1,1]=starting5new[1,4];grgamma00[1,2,2]=starting5new[2,4];grgamma00[1,1,3:11]=starting5new[3:11,4];grgamma00[1,2,12:20]=starting5new[12:20,4]
grgamma00[2,1,1]=starting5new[1,5];grgamma00[2,2,2]=starting5new[2,5];grgamma00[2,1,3:11]=starting5new[3:11,5];grgamma00[2,2,12:20]=starting5new[12:20,5]
grbeta00=as.matrix(starting5new[,6:7])

mu100=c(0,0)
mu200=c(0.1472920, -0.0157996)
mu300=c(0.1100240, -0.1232052)
Sig100=matrix(c(1,0.8512375,0.8512375,1),2,2)
Sig200=matrix(c(1,0.846,0.846,1),2,2)
Sig300=matrix(c(1,0.858,0.858,1),2,2)


# 1. EM lasso method

if (min(resp)==0){
  resp2=as.matrix(resp)
  resp=resp+1
} else {
  resp2=as.matrix(resp)-1
}
# tuning parameter
eta=35
ptm <- proc.time()
sim=Reg_DIF(resp=resp,m=2,r=2,y=3,N.vec=c(N1,N2,N3),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
print(proc.time() - ptm)

# 1. EM lasso method

# tuning parameter
eta=21
ptm <- proc.time()
sim=Reg_DIF(resp=resp,m=m,r=r,y=y,N.vec=c(N1,N2,N3),eta=eta,eps =1e-3,max.tol=1e-7,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
print(proc.time() - ptm)
est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta,iter=iter,bic=BIC, means=Mu.est,Covs=Sig.est
a1=sim$grgamma[c(4,5,12,13),1])!=0
#(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(4,5,12,13),1])!=0
a2=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(2,3,10,11),1])!=0
a3=(read.csv(paste("/Users/zhux0445/Documents/GitHub/RegDIF_SimData/NABeta5_",rep),row.names = 1)[c(4,5,12,13),2])!=0
a4=(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/NAGamma5_",rep),row.names = 1)[c(2,3,10,11),2])!=0

a5=(a1|a2)!=0
a6=(a3|a4)!=0

# 2. EMM lasso method

ptm <- proc.time()
sim2=Reg_EMM_DIF(resp,m,r,y,c(N1,N2,N3),eta,eps =1e-3,max.tol=1e-7,gra00,grd00,grbeta00,grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
print(proc.time() - ptm)
sim2

# 3. EM Adaptive lasso method
eta=13
ptm <- proc.time()
sim3=Reg_Adaptive_DIF(resp,m,r,y,c(N1,N2,N3),eta,lam=1,eps =1e-3,max.tol=1e-7,gra00,grd00,grbeta00,grgamma00,Mu.list=c(mu100,mu200,mu300),Sig.list=rbind(Sig100,Sig200,Sig300))
print(proc.time() - ptm)
sim3
  