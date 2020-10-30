library(mirt) 
library(MASS)
library(mvtnorm)
mirtCluster(4)
# uniform DIF, 2D2PL 
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para4.csv",row.names = 1)
Amat1=as.matrix(params[,1:2])
rownames(Amat1) <- c()
Amat2=as.matrix(params[,3:4])
rownames(Amat2) <- c()
Amat2=as.matrix(params[,5:6])
rownames(Amat3) <- c()
Dmat1=matrix(params[,7],J,(m-1))
Dmat2=matrix(params[,8],J,(m-1))
Dmat2=matrix(params[,9],J,(m-1))


# Non-uniform DIF on slope only, 2D2PL 


# calculating wABC to check generated DIF magnitude
X1=seq(-4,4,by=0.25)
r=2
G=length(X1)^r
gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
X <- matrix(gh[idx,1],nrow(idx),r)
WXr=WXf1=dmvnorm(X,c(0,0),matrix(c(1,0.85,0.85,1),2,2))
#WXf2=dmvnorm(X,c(0,0),matrix(c(1,0.15,0.15,1),2,2)) # for impact condition 
G=length(X1)^r
m=2
J=20
gra=Amat1
grd=Dmat1
ygam1=matrix(0,J,2)
ygam2=Amat2-Amat1
ygamat1=ygam1%*%t(X)
ygamat2=ygam2%*%t(X)
ybeta1=matrix(0,J,1)
ybeta2=Dmat2-Dmat1
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
  pstar2[,,g] = 1/(1+exp(-(grd+matrix(axmat[,g],J,1)+matrix(ygamat2[,g],J,1)+ybeta2))) #(3)
  p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
}
wABCr=wABCf1=wABCf2=numeric(J)
for (j in 1:J){
  wABCr[j]=sum(abs((p1[,2,]-p2[,2,])[j,])*WXr)
  wABCf1[j]=sum(abs((p1[,2,]-p2[,2,])[j,])*WXf1)
  #wABCf2[j]=sum(abs((p1[,2,]-p2[,2,])[j,])*WXf2)
}
wABC1=  (wABCr+  wABCf1)/2 #for no impact condition, item 6,11,12,13,14 all have wABC between 0.3-0.9
#wABC2=  (wABCr+  wABCf2)/2 #for with impact condition, item 6,11,12,13,14 all have wABC between 0.3-0.9

