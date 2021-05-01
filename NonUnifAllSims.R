#MBR non-uniform DIF simulations
library(mirt) 
library(MASS)
library(mvtnorm)
mirtCluster(4)
library(ggplot2)

# uniform DIF, 2D2PL 
#set.seed(1)
#A=runif(20,1.5,2.5)
A11=c(2.17,0,2.41,2.45,2.34,1.84,1.85,1.92,1.94,1.90,1.92,rep(0,9))
A12=A11;A13=A11
A12[4:9]=A12[4:9]-0.4;A13[4:9]=A13[4:9]-0.6
A21=c(0,2.46,rep(0,9),2.43,1.82,2.22,1.93,1.88,1.84,2.12,2.42,2.15)
A22=A21;A23=A21
A22[12:17]=A22[12:17]-0.4;A23[12:17]=A23[12:17]-0.6
D11=c(0.03,-1.28,0.58,-2.06,0.12,3.25,-0.41,-0.51,0.89,1.33,0.85,0.82,-0.37,-0.99,-0.27,0.19,1.73,0.05,-1.86,-0.63)
D12=D11;D13=D11
D12[c(4:9,12:17)]=D12[c(4:9,12:17)]+0.25
D13[c(4:9,12:17)]=D13[c(4:9,12:17)]+0.6

Amat1=cbind(A11,A21);Amat2=cbind(A12,A22);Amat3=cbind(A13,A23)
Dmat1=cbind(D11)
Dmat2=cbind(D12);Dmat3=cbind(D13)


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
