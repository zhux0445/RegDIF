#BIC in Tutz

setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
# Dataset #4 (2 Non-uniform DIF items per scale)
J=20

set.seed(1)
A11=c(runif(1,1.5,2.5),0,runif(9,1.5,2.5),rep(0,9))
A12=A13=A11
A21=c(0,runif(1,1.5,2.5),rep(0,9),runif(9,1.5,2.5))
A22=A23=A21
D11=rnorm(J,0,1)
D12=D13=D11
# focal group has smaller slopes and more extreme intercepts
A12[4]=A11[4]- 0.3
A12[5]=A11[5]- 0.3
A22[12]=A21[12]- 0.3
A22[13]=A21[13]- 0.3
A13[4]=A11[4]- 0.6
A13[5]=A11[5]- 0.6
A23[12]=A21[12]- 0.6
A23[13]=A21[13]- 0.6
D12[4]=D11[4]+ 0.5
D12[5]=D11[5]+ 0.5
D13[4]=D11[4]+ 1
D13[5]=D11[5]+ 1
D12[12]=D11[12]+ 0.5
D12[13]=D11[13]+ 0.5
D13[12]=D11[12]+ 1
D13[13]=D11[13]+ 1
Amat1=cbind(A11,A21)
Amat2=cbind(A12,A22)
Amat3=cbind(A13,A23)
Dmat1=cbind(D11)
Dmat2=cbind(D12)
Dmat3=cbind(D13)

N1=N2=N3=500
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
N=N1+N2+N3

m=2
r=2
y1=c(0,0)
y2=c(1,0)
y3=c(0,1)
###########################################################
######################             ########################
######################    Seed     ########################
######################             ########################
###########################################################
set.seed(200)
Theta=mvrnorm(n=N,mu=c(0,0),Sigma = matrix(c(1,0.85,0.85,1),2,2))
Theta1=Theta[1:N1,]
Theta2=Theta[(N1+1):(N1+N2),]
Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
resp=rbind(datasetR,datasetF1,datasetF2)

########################################################
###                                                 ####
###    Group Lasso DIF Quadrature (mirt start)      ####
###                                                 ####
########################################################

library(mirt) 
mirtCluster(4)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(irtoys)
library(microbenchmark)
library(Rcpp)
sourceCpp("/Users/ruoyizhu/Documents/GitHub/mirt/matrix.cpp")

soft=function(s, tau) {      
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

########################################################
#    starting values by multiplegroup for all etas     #
########################################################
s <- 'F1 = 1,3-11
          F2 = 2,12-20
          COV = F1*F2'
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
md00 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:r]))
gra000=coef(md00,simplify=T)$G1$items[,c("a1","a2")]
grd000=matrix(coef(md00,simplify=T)$G1$items[,c("d")],20,1)
grgamma000=array(0,dim=c(r,r,J))
grgamma000[1,1,3:11]=(coef(md00,simplify=T)$G2$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[3:11,1]
grgamma000[2,1,3:11]=(coef(md00,simplify=T)$G3$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[3:11,1]
grgamma000[1,2,12:20]=(coef(md00,simplify=T)$G2$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[12:20,2]
grgamma000[2,2,12:20]=(coef(md00,simplify=T)$G3$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[12:20,2]
grbeta000=matrix(0,J,2)
grbeta000[,1]=matrix(coef(md00,simplify=T)$G2$items[,c("d")],20,1)-matrix(coef(md00,simplify=T)$G1$items[,c("d")],20,1)
grbeta000[,2]=matrix(coef(md00,simplify=T)$G3$items[,c("d")],20,1)-matrix(coef(md00,simplify=T)$G1$items[,c("d")],20,1)
Sig100=coef(md00,simplify=T)$G1$cov
Sig200=coef(md00,simplify=T)$G2$cov
Sig300=coef(md00,simplify=T)$G3$cov
Mu100=coef(md00,simplify=T)$G1$means
Mu200=coef(md00,simplify=T)$G2$means
Mu300=coef(md00,simplify=T)$G3$means


##############################################
#    starting values by mirt for eta=0       #
##############################################
s <- 'F1 = 1,3-11
          F2 = 2,12-20
          COV = F1*F2'
cmodel <- mirt.model(s)
md0 <- mirt(resp, cmodel, itemtype = "2PL")
gra0=coef(md0,simplify=T)$items[,1:r]
grd0=matrix(coef(md0,simplify=T)$items[,3],J,1)
grgamma0=array(0,dim=c(r,r,J))
grgamma0[1,1,3:11]=-0.7
grgamma0[2,1,3:11]=-1.2
grgamma0[1,2,12:20]=-0.7
grgamma0[2,2,12:20]=-1.2
grbeta0=matrix(0,J,2)
grbeta0[,1]=matrix(c(0,0,rep(0.7,18)),J,1)
grbeta0[,2]=matrix(c(0,0,rep(1.4,18)),J,1)
Sig0=coef(md0,simplify=T)$cov

##############################################
#                                            #
#  starting values at eta=0 for later eta's  #
#                                            #
##############################################
eta=0
eps =1e-3
max.tol=1e-7
NonUniform=T
if(min(resp)==0)
{ 
  resp <- resp+1
}
#sample size and test length
N <- nrow(resp)
J <- ncol(resp)

X1=seq(-3,3,by=0.2)
G=length(X1)^r
gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
X <- matrix(gh[idx,1],nrow(idx),r)
ng = ng1 = ng2 = ng3 = numeric(G)
Xijk=array(double(N*J*m),dim = c(N,J,m))
for(i in 1:N){
  for(j in 1:J){
    for(k in 1:m){
      Xijk[i,j,k]=ifelse(resp[i,j]==k,1,0)
    }
  }
}
#gra = azero
#grd = matrix(0,J,1)
gra=gra0
grd=grd0
#grbeta=grbeta0
Sig.gp1=Sig.gp2=Sig.gp3=Sig0
grgamma=grgamma0
grbeta=grbeta0
#starting values for the MC samples form the prior
Mu.gp1=Mu.gp2=Mu.gp3=rep(0,r)
#X=theta #for check later

Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
df.a <- df.d  <- df.gamma  <- df.beta  <- 1

iter <- 0

while(max(df.a)>eps | max(df.d)>eps |  max(df.gamma)>eps | max(df.beta)>eps)
{
  aold <- gra
  dold <- grd
  #gammaold=grgamma
  gammaold=grgamma
  betaold=grbeta
  A1=dmvnorm(X,Mu.gp1,Sig.gp1)
  A2=dmvnorm(X,Mu.gp2,Sig.gp2)
  A3=dmvnorm(X,Mu.gp3,Sig.gp3)
  
  #calculation of n_g 
  axmat=eigenMapMatMult(gra,t(X))
  ygam1=ygam2=ygam3=matrix(0,J,r)
  ygam2=t(grgamma[1,,])
  ygam3=t(grgamma[2,,])
  ygamat1=eigenMapMatMult(ygam1,t(X))
  ygamat2=eigenMapMatMult(ygam2,t(X))
  ygamat3=eigenMapMatMult(ygam3,t(X))
  grbeta1=t(y1%*%t(grbeta))
  grbeta2=t(y2%*%t(grbeta))
  grbeta3=t(y3%*%t(grbeta))
  pstar1=pstar2=pstar3=array(double(J*(m-1)*G),dim = c(J,(m-1),G))
  p1=p2=p3=array(double(J*m*G),dim = c(J,m,G))
  for (g in 1:G)
  {
    pstar1[,,g] = 1/(1+exp(-(grd+axmat[,g]+ygamat1[,g]+grbeta1)))
    p1[,,g] = t(-diff(rbind(rep(1,J),t(pstar1[,,g]),rep(0,J))))
  }
  for (g in 1:G)
  {
    pstar2[,,g] = 1/(1+exp(-(grd+axmat[,g]+ygamat2[,g]+grbeta2)))
    p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
  }
  for (g in 1:G)
  {
    pstar3[,,g] = 1/(1+exp(-(grd+axmat[,g]+ygamat3[,g]+grbeta3)))
    p3[,,g] = t(-diff(rbind(rep(1,J),t(pstar3[,,g]),rep(0,J))))
  }
  
  pij=array(double(J*G*N1),dim=c(J,G,N1))
  LiA=matrix(double(N*G),N,G)
  
  for (j in 1:J)
  {
    for (g in 1:G)
    {
      pij[j,g,]=p1[j,(resp[1:N1,j]),g]
    }
  }
  LiA[1:N1,]=t(apply(pij, c(2,3), prod)*A1)
  
  for (j in 1:J)
  {
    for (g in 1:G)
    {
      pij[j,g,]=p2[j,(resp[(N1+1):(N1+N2),j]),g]
    }
  }
  LiA[(N1+1):(N1+N2),]=t(apply(pij, c(2,3), prod)*A2)
  
  for (j in 1:J)
  {
    for (g in 1:G)
    {
      pij[j,g,]=p3[j,(resp[(N1+N2+1):(N1+N2+N3),j]),g]
    }
  }
  LiA[(N1+N2+1):(N1+N2+N3),]=t(apply(pij, c(2,3), prod)*A3)
  
  
  
  Pi = apply(LiA,1,sum)
  ng = apply(LiA/Pi,2,sum)
  ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
  ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
  ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
  
  #update mu hat and Sigma hat
  Mu.gp2=colSums(X*ng2)/N2
  Mu.gp3=colSums(X*ng3)/N3
  
  #update mu hat and Sigma hat
  Sig.hat1=eigenMapMatMult(t(X),(X*ng1))/N1
  Sig.hat2=eigenMapMatMult(t(X-rep(Mu.gp2)),((X-rep(Mu.gp2))*ng2))/N2
  Sig.hat3=eigenMapMatMult(t(X-rep(Mu.gp3)),((X-rep(Mu.gp3))*ng3))/N3
  
  #scale 
  #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
  Tau=sqrt(diag(Sig.hat1))
  Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
  #q_g_star
  #X=(X-mu.hat.mat)/Tau.mat
  Xstar=X/Tau.mat
  
  Sig.gp1=eigenMapMatMult(t(Xstar),(Xstar*ng1))/N1
  Sig.gp2=eigenMapMatMult(t(Xstar-rep(Mu.gp2)),((Xstar-rep(Mu.gp2))*ng2))/N2
  Sig.gp3=eigenMapMatMult(t(Xstar-rep(Mu.gp3)),((Xstar-rep(Mu.gp3))*ng3))/N3
  
  
  ##Constraint part
  #calculation of r_jgk
  for (j in 1:r)
  {
    d <- grd[j,]
    a <- gra[j,j]
    rLiA <- array(double(N*G*m),dim = c(N,G,m))
    for(k in 1:m){
      rLiA[,,k]=Xijk[,j,k]*LiA/Pi
    }
    rgk <- apply(rLiA,c(2,3),sum)
    
    #M-step loop starts for item j
    miter <- 0
    add <- max.tol+1
    while(sum(abs(add))>max.tol)
    {
      miter <- miter+1
      for(g in 1:G){
        Pstar[g,] <- 1/(1+exp(-(d+rep(a*X[g,j]))))
        P[g,] <- -diff(c(1,Pstar[g,],0))
      }
      Qstar <- 1-Pstar
      
      #calculating the score vector
      if (m==2){
        Dsco <- sum(apply(rgk/P,1,diff)*Pstar*Qstar)
      }else
        Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      
      PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
      Asco <- sum(X[,j]*apply(rgk*PQdif/P,1,sum))
      minusgrad <- -c(Dsco,Asco)
      
      FI <- matrix(0,m,m)
      if (m==2){
        for (mm in 1:(m-1)){
          FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
          FI[m,mm] <- FI[mm,m] <- sum(ng*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
        }
        
      } else {
        for (mm in 1:(m-1)){
          FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
          FI[m,mm] <- FI[mm,m] <- sum(ng*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
        }
        for (mm in 2:(m-1)){
          FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
        }
        for (mm in 1:(m-2)){
          FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
        }
      }
      
      FI[m,m] <- -sum(ng*X[,j]^2*apply(PQdif^2/P,1,sum))
      
      
      add <- qr.solve(FI,minusgrad)
      
      d <- d+add[1:(m-1)]
      a <- a+add[m]
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[j]
    gra[j,j] <- a
    grd[j,] <- d
  }   
  
  
  ##random part
  # dimension1
  #calculation of r_jgk
  for (j in (r+1):(r+J/2-1))
  {
    d <- grd[j,] 
    a <- gra[j,1]
    gam=grgamma[,,j]
    bet=grbeta[j,]
    rLiA <- array(double(N*G*m),dim = c(N,G,m))
    for(k in 1:m){
      rLiA[,,k]=Xijk[,j,k]*LiA/Pi
    }
    rgk <- apply(rLiA,c(2,3),sum)
    rgk1 <- apply(rLiA[1:N1,,],c(2,3),sum)
    rgk2 <- apply(rLiA[(N1+1):(N1+N2),,],c(2,3),sum)
    rgk3 <- apply(rLiA[(N1+N2+1):(N1+N2+N3),,],c(2,3),sum)
    
    #M-step loop starts for item j
    miter <- 0
    add <- max.tol+1
    while(sum(abs(add))>max.tol)
    {
      miter <- miter+1
      for(g in 1:G){
        Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y1%*%bet+rep((y1%*%gam)%*%X[g,]))))
        P1[g,] <- -diff(c(1,Pstar1[g,],0))
        Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y2%*%bet+rep((y2%*%gam)%*%X[g,]))))
        P2[g,] <- -diff(c(1,Pstar2[g,],0))
        Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y3%*%bet+rep((y3%*%gam)%*%X[g,]))))
        P3[g,] <- -diff(c(1,Pstar3[g,],0))
      }
      Qstar <- 1-Pstar
      Qstar1 <- 1-Pstar1
      Qstar2 <- 1-Pstar2
      Qstar3 <- 1-Pstar3
      
      
      #calculating the score vector
      if (m==2){
        Dsco <- sum(apply(rgk1/P1,1,diff)*Pstar1*Qstar1+apply(rgk2/P2,1,diff)*Pstar2*Qstar2+apply(rgk3/P3,1,diff)*Pstar3*Qstar3)
      }else
        Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
      PQdif1 <- -t(apply(cbind(0,Pstar1*Qstar1,0),1, diff))
      PQdif2 <- -t(apply(cbind(0,Pstar2*Qstar2,0),1, diff))
      PQdif3 <- -t(apply(cbind(0,Pstar3*Qstar3,0),1, diff))
      Asco = sum(X[,1]*apply(rgk1/P1*PQdif1,1,sum))+sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum))+sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum))
      if (NonUniform==T){
        Gamsco = c(sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum)))
      }
      if (m==2){
        Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      }else
        Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      if (NonUniform==T){
        minusgrad <- -c(Dsco,Asco,Gamsco,Betsco)# ,Betsco
        grad <- c(Dsco,Asco,Gamsco,Betsco)# ,Betsco
      } else {
        minusgrad <- -c(Dsco,Asco, Betsco)
        grad <- c(Dsco,Asco, Betsco)
      }
      
      if (NonUniform==T){
        FI <- matrix(0,length(minusgrad),length(minusgrad))#FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
          }
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
            }
          }
          if (m>2){
            for (mm in 2:(m-1)){
              FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
            }
            for (mm in 1:(m-2)){
              FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
            }
          }
        }
        FI[m,m] <- -sum(ng1*X[,1]*X[,1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
        FI[3,3]=FI[3,2]=FI[2,3]=-sum(ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum))
        FI[4,4]=FI[4,2]=FI[2,4]=-sum(ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
        FI[3,1]=FI[1,3]=sum(ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1]))
        FI[4,1]=FI[1,4]=sum(ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
        for (mm in 5:6){
          #(5,5) (6,6)
          FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
          #(2,5) (5,2) (2,6) (6,2)
          FI[m,mm] <- FI[mm,m] <- FI[mm,mm-2] <- FI[mm-2,mm] <- c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-4]
          #(1,5) (5,1) (1,6) (6,1)
          FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
        }
      } else {
        FI <- matrix(0,m+(m-1)*2,m+(m-1)*2)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
          }
          for (mm in (m+1):(m+(m-1)*2)){
            FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
            FI[m,mm] <- FI[mm,m] <- c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-m]
            FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
          }
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
            }
          }
          if (m>2){
            for (mm in 2:(m-1)){
              FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
            }
            for (mm in 1:(m-2)){
              FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
            }
          }
        }
        FI[m,m] <- -sum(ng1*X[,1]*X[,1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
      }
      #sqrt(sum((grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6])^2))
      #gp=grad-c(d,a,gam[,1],bet)%*%FI
      gp=grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6]
      l2gp=sqrt(innerProduct(grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6],grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6]))
      if (l2gp<eta){
        gam[,1]=bet[1:2]=0
        add <- qr.solve(FI,minusgrad)
        add[3:6]=0
        d=d+add[1:(m-1)]
        a=a+add[m]
      } else {
        add <- qr.solve(FI,minusgrad)
        direction=qr.solve(FI[3:6,3:6],as.vector(eta*gp/l2gp-grad[3:6]))
        #for (ii in 1:4){
        #  if (abs(direction[ii])<1e-3){
        #    direction[ii]=0
        #  }
        #}
        #direction=zapsmall(c(c(gam[,1],bet),direction), digits = 9)[5:8]
        linesch.alpha=1
        #if (direction[1] !=0 |direction[2] !=0 |direction[3] !=0 |direction[4] !=0 )
        #Armijo line search
        #{
        #deltaj= direction%*%minusgrad[3:6]+eta*sqrt(innerProduct(c(gam[,1],bet)+direction,c(gam[,1],bet)+direction))-eta*sqrt(innerProduct(c(gam[,1],bet),c(gam[,1],bet)))
        #steplength.vec=c(0.5^c(0:10))
        #linesch=matrix(0,length(steplength.vec),4)
        #for (ii in 1:length(steplength.vec)){
        #  linesch[ii,]=direction*steplength.vec[ii]
        #}
        #updtedtau.vec=matrix(rep(c(gam[,1],bet),length(steplength.vec)),length(steplength.vec),4,byrow=T)+linesch
        #likelihood.vec=numeric(length(steplength.vec))
        #for (ii in 1:length(steplength.vec)){
        #  sumoverk1=numeric(G)
        #  for(g in 1:G){
        #    sumoverk1[g]=rgk1[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1])))),0)))
        #  }
        #  sumoverk2=numeric(G)
        #  for(g in 1:G){
        #    sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(updtedtau.vec[ii,1]%*%X[g,1])+ updtedtau.vec[ii,3]))),0)))
        #  }
        #  sumoverk3=numeric(G)
        #  for(g in 1:G){
        #    sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(updtedtau.vec[ii,2]%*%X[g,1])+ updtedtau.vec[ii,4]))),0)))
        #  }
        #  likelihood.vec[ii]=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)
        #}
        #sumoverk2=numeric(G)
        #for(g in 1:G){
        #  sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(gam[1,1]%*%X[g,1])+ bet[1]))),0)))
        #}
        #sumoverk3=numeric(G)
        #for(g in 1:G){
        #  sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(gam[2,1]%*%X[g,1])+ bet[2]))),0)))
        #}
        #lh0=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)
        #Armijols=numeric(length(steplength.vec))
        #for (ii in 1:length(steplength.vec)){
        #  Armijols[ii]=(lh0-likelihood.vec[ii]+eta*sqrt(innerProduct(c(gam[,1],bet)+direction*steplength.vec[ii],c(gam[,1],bet)+direction*steplength.vec[ii]))-eta*sqrt(innerProduct(c(gam[,1],bet),c(gam[,1],bet)))<=steplength.vec[ii]*0.1* deltaj)
        #}
        #if (max(Armijols)==0){
        #  linesch.alpha=0
        #} else {
        #  linesch.alpha=max(steplength.vec[which(Armijols==1)])
        #  }
        #}
        d=d+add[1:(m-1)]
        a=a+add[m]
        add[3:6]=linesch.alpha*direction
        gam[,1]=gam[,1]+add[3:4]
        bet=bet+add[5:6]
      }
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[1]
    gam=gam*Tau[1]
    
    gra[j,1] <- a
    grd[j,] <- d
    grgamma[,,j] <-gam
    grbeta[j,]=bet
  }
  
  # dimension2
  for (j in (r+J/2):J)
  {
    d <- grd[j,] 
    a <- gra[j,2]
    gam=grgamma[,,j]
    bet=grbeta[j,]
    rLiA <- array(double(N*G*m),dim = c(N,G,m))
    for(k in 1:m){
      rLiA[,,k]=Xijk[,j,k]*LiA/Pi
    }
    rgk <- apply(rLiA,c(2,3),sum)
    rgk1 <- apply(rLiA[1:N1,,],c(2,3),sum)
    rgk2 <- apply(rLiA[(N1+1):(N1+N2),,],c(2,3),sum)
    rgk3 <- apply(rLiA[(N1+N2+1):(N1+N2+N3),,],c(2,3),sum)
    
    #M-step loop starts for item j
    miter <- 0
    add <- max.tol+1
    while(sum(abs(add))>max.tol)
    {
      miter <- miter+1
      for(g in 1:G){
        Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y1%*%bet+rep((y1%*%gam)%*%X[g,]))))
        P1[g,] <- -diff(c(1,Pstar1[g,],0))
        Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y2%*%bet+rep((y2%*%gam)%*%X[g,]))))
        P2[g,] <- -diff(c(1,Pstar2[g,],0))
        Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y3%*%bet+rep((y3%*%gam)%*%X[g,]))))
        P3[g,] <- -diff(c(1,Pstar3[g,],0))
      }
      Qstar <- 1-Pstar
      Qstar1 <- 1-Pstar1
      Qstar2 <- 1-Pstar2
      Qstar3 <- 1-Pstar3
      
      
      #calculating the score vector
      if (m==2){
        Dsco <- sum(apply(rgk1/P1,1,diff)*Pstar1*Qstar1+apply(rgk2/P2,1,diff)*Pstar2*Qstar2+apply(rgk3/P3,1,diff)*Pstar3*Qstar3)
      }else
        Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
      PQdif1 <- -t(apply(cbind(0,Pstar1*Qstar1,0),1, diff))
      PQdif2 <- -t(apply(cbind(0,Pstar2*Qstar2,0),1, diff))
      PQdif3 <- -t(apply(cbind(0,Pstar3*Qstar3,0),1, diff))
      Asco = sum(X[,2]*apply(rgk1/P1*PQdif1,1,sum))+sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum))+sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum))
      if (NonUniform==T){
        Gamsco = c(sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum)))
      }
      if (m==2){
        Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      }else
        Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      if (NonUniform==T){
        minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)#, Betsco
        grad <- c(Dsco,Asco,Gamsco,Betsco)#, Betsco
      } else {
        minusgrad <- -c(Dsco,Asco, Betsco)
        grad <- c(Dsco,Asco, Betsco)
      }
      
      if (NonUniform==T){
        FI <- matrix(0,length(minusgrad),length(minusgrad))#FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,2]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
          }
        } else { 
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
            }
          }
          if (m>2){
            for (mm in 2:(m-1)){
              FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
            }
            for (mm in 1:(m-2)){
              FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
            }
          }
        }
        FI[m,m] <- -sum(ng1*X[,2]*X[,2]*apply(PQdif1^2/P1,1,sum)+ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum)+ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
        FI[3,3]=FI[3,2]=FI[2,3]=-sum(ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum))
        FI[4,4]=FI[4,2]=FI[2,4]=-sum(ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
        FI[3,1]=FI[1,3]=sum(ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1]))
        FI[4,1]=FI[1,4]=sum(ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
        for (mm in 5:6){
          #(5,5) (6,6)
          FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
          #(2,5) (5,2) (2,6) (6,2)
          FI[m,mm] <- FI[mm,m] <- FI[mm,mm-2] <- FI[mm-2,mm] <- c(sum(ng2*X[,2]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,2]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-4]
          #(1,5) (5,1) (1,6) (6,1)
          FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
        }
      } else {
        FI <- matrix(0,m+(m-1)*2,m+(m-1)*2)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,2]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
          }
          for (mm in (m+1):(m+(m-1)*2)){
            FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
            FI[m,mm] <- FI[mm,m] <- c(sum(ng2*X[,2]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,2]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-m]
            FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
          }
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
            }
          }
          if (m>2){
            for (mm in 2:(m-1)){
              FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
            }
            for (mm in 1:(m-2)){
              FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
            }
          }
        }
        FI[m,m] <- -sum(ng1*X[,2]*X[,2]*apply(PQdif1^2/P1,1,sum)+ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum)+ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
      }
      
      gp=grad[3:6]-c(gam[,2],bet)%*%FI[3:6,3:6]
      l2gp=sqrt(innerProduct(grad[3:6]-c(gam[,2],bet)%*%FI[3:6,3:6],grad[3:6]-c(gam[,2],bet)%*%FI[3:6,3:6]))
      if (l2gp<eta){
        gam[,2]=bet[1:2]=0
        add <- qr.solve(FI,minusgrad)
        add[3:6]=0
        d=d+add[1:(m-1)]
        a=a+add[m]
      } else {
        add <- qr.solve(FI,minusgrad)
        direction=qr.solve(FI[3:6,3:6],as.vector(eta*gp/l2gp-grad[3:6]))
        #for (ii in 1:4){
        #  if (abs(direction[ii])<1e-3){
        #    direction[ii]=0
        #  }
        #}
        #direction=zapsmall(c(c(gam[,2],bet),direction), digits = 9)[5:8]
        linesch.alpha=1
        #if (direction[1] !=0 |direction[2] !=0 |direction[3] !=0 |direction[4] !=0 )
        #{
        #  deltaj= direction%*%minusgrad[3:6]+eta*sqrt(innerProduct(c(gam[,2],bet)+direction,c(gam[,2],bet)+direction))-eta*sqrt(innerProduct(c(gam[,2],bet),c(gam[,2],bet)))
        #  steplength.vec=c(0.5^c(0:10))
        #  linesch=matrix(0,length(steplength.vec),4)
        #  for (ii in 1:length(steplength.vec)){
        #    linesch[ii,]=direction*steplength.vec[ii]
        #  }
        #  updtedtau.vec=matrix(rep(c(gam[,2],bet),length(steplength.vec)),length(steplength.vec),4,byrow=T)+linesch
        #  likelihood.vec=numeric(length(steplength.vec))
        #  for (ii in 1:length(steplength.vec)){
        #    sumoverk2=numeric(G)
        #    for(g in 1:G){
        #      sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(updtedtau.vec[ii,1]%*%X[g,2])+ updtedtau.vec[ii,3]))),0)))
        #    }
        #    sumoverk3=numeric(G)
        #    for(g in 1:G){
        #      sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(updtedtau.vec[ii,2]%*%X[g,2])+ updtedtau.vec[ii,4]))),0)))
        #    }
        #    likelihood.vec[ii]=sum(sumoverk2)+sum(sumoverk3)
        #  }
        #  sumoverk2=numeric(G)
        #  for(g in 1:G){
        #    sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(gam[1,2]%*%X[g,2])+ bet[1]))),0)))
        #  }
        #  sumoverk3=numeric(G)
        #  for(g in 1:G){
        #    sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(gam[2,2]%*%X[g,2])+ bet[2]))),0)))
        #  }
        #  lh0=sum(sumoverk2)+sum(sumoverk3)
        #  Armijols=numeric(length(steplength.vec))
        #  for (ii in 1:length(steplength.vec)){
        #    Armijols[ii]=(lh0-likelihood.vec[ii]+eta*sqrt(innerProduct(c(gam[,2],bet)+direction*steplength.vec[ii],c(gam[,2],bet)+direction*steplength.vec[ii]))-eta*sqrt(innerProduct(c(gam[,2],bet),c(gam[,2],bet)))<=steplength.vec[ii]*0.1* deltaj)
        #  }
        #  if (max(Armijols)==0){
        #    linesch.alpha=0
        #  } else {
        #    linesch.alpha=max(steplength.vec[which(Armijols==1)])
        #  }
        #}
        d=d+add[1:(m-1)]
        a=a+add[m]
        add[3:6]=linesch.alpha*direction
        #add[3:6]=direction
        gam[,2]=gam[,2]+add[3:4]
        bet=bet+add[5:6]
      }
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[2]
    gam=gam*Tau[2]
    
    gra[j,2] <- a
    grd[j,] <- d
    grgamma[,,j] <-gam
    grbeta[j,]=bet
  }
  
  df.d <- abs(dold-grd)
  df.a <- abs(aold-gra)
  df.gamma <- abs(gammaold-grgamma)
  df.beta = abs(betaold-grbeta)
  iter <- iter+1
}


gra00=gra
grd00=grd
grgamma00=grgamma
grbeta00=grbeta
mu100= Mu.gp1
mu200= Mu.gp2
mu300= Mu.gp3
Sig100=Sig.gp1
Sig200=Sig.gp2
Sig300=Sig.gp3

#########################
##     For \eta>0      ##
#########################

########################
####### r=2, m=2 #######
########################

#########################
##                     ##
##  start of function  ##
##                     ##
#########################

ipest1 <- function(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=T,gra00=gra000,grd00=grd000,grgamma00=grgamma000,grbeta00=grbeta000,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
{
  #make sure the responses are coded from 1 instead of 0
  if(min(resp)==0)
  {
    resp <- resp+1
  }
  #sample size and test length
  N <- nrow(resp)
  J <- ncol(resp)
  
  X1=seq(-3,3,by=0.2)
  G=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X <- matrix(gh[idx,1],nrow(idx),r)
  ng= ng1= ng2= ng3 = numeric(G)
  Xijk=array(double(N*J*m),dim = c(N,J,m))
  for(i in 1:N){
    for(j in 1:J){
      for(k in 1:m){
        Xijk[i,j,k]=ifelse(resp[i,j]==k,1,0)
      }
    }
  }
  #gra = azero
  #grd = matrix(0,J,1)
  gra=gra00
  grd=grd00
  grgamma=grgamma00
  grbeta=grbeta00
  Sig.gp1=Sig100; Sig.gp2=Sig200; Sig.gp3=Sig300
  
  #starting values for the MC samples form the prior
  Mu.gp1=mu100; Mu.gp2=mu200; Mu.gp3=mu300
  #X=theta #for check later
  
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.gamma <- df.beta <- 1
  
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps | max(df.d)>eps |  max(df.gamma)>eps | max(df.beta)>eps)
  {
    aold <- gra
    dold <- grd
    #gammaold=grgamma
    gammaold=grgamma
    betaold=grbeta
    A1=dmvnorm(X,Mu.gp1,Sig.gp1)
    A2=dmvnorm(X,Mu.gp2,Sig.gp2)
    A3=dmvnorm(X,Mu.gp3,Sig.gp3)
    
    #calculation of n_g 
    axmat=eigenMapMatMult(gra,t(X))
    ygam1=ygam2=ygam3=matrix(0,J,r)
    ygam2=t(grgamma[1,,])
    ygam3=t(grgamma[2,,])
    ygamat1=eigenMapMatMult(ygam1,t(X))
    ygamat2=eigenMapMatMult(ygam2,t(X))
    ygamat3=eigenMapMatMult(ygam3,t(X))
    grbeta1=t(y1%*%t(grbeta))
    grbeta2=t(y2%*%t(grbeta))
    grbeta3=t(y3%*%t(grbeta))
    pstar1=pstar2=pstar3=array(double(J*(m-1)*G),dim = c(J,(m-1),G))
    p1=p2=p3=array(double(J*m*G),dim = c(J,m,G))
    for (g in 1:G)
    {
      pstar1[,,g] = 1/(1+exp(-(grd+axmat[,g]+ygamat1[,g]+grbeta1)))
      p1[,,g] = t(-diff(rbind(rep(1,J),t(pstar1[,,g]),rep(0,J))))
    }
    for (g in 1:G)
    {
      pstar2[,,g] = 1/(1+exp(-(grd+axmat[,g]+ygamat2[,g]+grbeta2)))
      p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
    }
    for (g in 1:G)
    {
      pstar3[,,g] = 1/(1+exp(-(grd+axmat[,g]+ygamat3[,g]+grbeta3)))
      p3[,,g] = t(-diff(rbind(rep(1,J),t(pstar3[,,g]),rep(0,J))))
    }
    
    pij=array(double(J*G*N1),dim=c(J,G,N1))
    LiA=matrix(double(N*G),N,G)
    
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g,]=p1[j,(resp[1:N1,j]),g]
      }
    }
    LiA[1:N1,]=t(apply(pij, c(2,3), prod)*A1)
    
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g,]=p2[j,(resp[(N1+1):(N1+N2),j]),g]
      }
    }
    LiA[(N1+1):(N1+N2),]=t(apply(pij, c(2,3), prod)*A2)
    
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g,]=p3[j,(resp[(N1+N2+1):(N1+N2+N3),j]),g]
      }
    }
    LiA[(N1+N2+1):(N1+N2+N3),]=t(apply(pij, c(2,3), prod)*A3)
    
    
    
    Pi = apply(LiA,1,sum)
    ng = apply(LiA/Pi,2,sum)
    ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
    ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
    ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
    
    #update mu hat and Sigma hat
    Mu.gp2=colSums(X*ng2)/N2
    Mu.gp3=colSums(X*ng3)/N3
    
    #update mu hat and Sigma hat
    Sig.hat1=eigenMapMatMult(t(X),(X*ng1))/N1
    Sig.hat2=eigenMapMatMult(t(X-rep(Mu.gp2)),((X-rep(Mu.gp2))*ng2))/N2
    Sig.hat3=eigenMapMatMult(t(X-rep(Mu.gp3)),((X-rep(Mu.gp3))*ng3))/N3
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat1))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    
    Sig.gp1=eigenMapMatMult(t(Xstar),(Xstar*ng1))/N1
    Sig.gp2=eigenMapMatMult(t(Xstar-rep(Mu.gp2)),((Xstar-rep(Mu.gp2))*ng2))/N2
    Sig.gp3=eigenMapMatMult(t(Xstar-rep(Mu.gp3)),((Xstar-rep(Mu.gp3))*ng3))/N3
    
    
    ##Constraint part
    #calculation of r_jgk
    for (j in 1:r)
    {
      d <- grd[j,]
      a <- gra[j,j]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(k in 1:m){
        rLiA[,,k]=Xijk[,j,k]*LiA/Pi
      }
      rgk <- apply(rLiA,c(2,3),sum)
      
      #M-step loop starts for item j
      miter <- 0
      add <- max.tol+1
      while(sum(abs(add))>max.tol)
      {
        miter <- miter+1
        for(g in 1:G){
          Pstar[g,] <- 1/(1+exp(-(d+rep(a*X[g,j]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        #calculating the score vector
        if (m==2){
          Dsco <- sum(apply(rgk/P,1,diff)*Pstar*Qstar)
        }else
          Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco <- sum(X[,j]*apply(rgk*PQdif/P,1,sum))
        minusgrad <- -c(Dsco,Asco)
        
        FI <- matrix(0,m,m)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          for (mm in 2:(m-1)){
            FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
          }
          for (mm in 1:(m-2)){
            FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
          }
        }
        
        FI[m,m] <- -sum(ng*X[,j]^2*apply(PQdif^2/P,1,sum))
        
        
        add <- qr.solve(FI,minusgrad)
        
        d <- d+add[1:(m-1)]
        a <- a+add[m]
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[j]
      gra[j,j] <- a
      grd[j,] <- d
    }   
    
    
    ##random part
    # dimension1
    #calculation of r_jgk
    for (j in (r+1):(r+J/2-1))
    {
      d <- grd[j,] 
      a <- gra[j,1]
      gam=grgamma[,,j]
      bet=grbeta[j,]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(k in 1:m){
        rLiA[,,k]=Xijk[,j,k]*LiA/Pi
      }
      rgk <- apply(rLiA,c(2,3),sum)
      rgk1 <- apply(rLiA[1:N1,,],c(2,3),sum)
      rgk2 <- apply(rLiA[(N1+1):(N1+N2),,],c(2,3),sum)
      rgk3 <- apply(rLiA[(N1+N2+1):(N1+N2+N3),,],c(2,3),sum)
      
      #M-step loop starts for item j
      miter <- 0
      add <- max.tol+1
      while(sum(abs(add))>max.tol)
      {
        miter <- miter+1
        for(g in 1:G){
          Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y1%*%bet+rep((y1%*%gam)%*%X[g,]))))
          P1[g,] <- -diff(c(1,Pstar1[g,],0))
          Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y2%*%bet+rep((y2%*%gam)%*%X[g,]))))
          P2[g,] <- -diff(c(1,Pstar2[g,],0))
          Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y3%*%bet+rep((y3%*%gam)%*%X[g,]))))
          P3[g,] <- -diff(c(1,Pstar3[g,],0))
        }
        Qstar <- 1-Pstar
        Qstar1 <- 1-Pstar1
        Qstar2 <- 1-Pstar2
        Qstar3 <- 1-Pstar3
        
        
        #calculating the score vector
        if (m==2){
          Dsco <- sum(apply(rgk1/P1,1,diff)*Pstar1*Qstar1+apply(rgk2/P2,1,diff)*Pstar2*Qstar2+apply(rgk3/P3,1,diff)*Pstar3*Qstar3)
        }else
          Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        PQdif1 <- -t(apply(cbind(0,Pstar1*Qstar1,0),1, diff))
        PQdif2 <- -t(apply(cbind(0,Pstar2*Qstar2,0),1, diff))
        PQdif3 <- -t(apply(cbind(0,Pstar3*Qstar3,0),1, diff))
        Asco = sum(X[,1]*apply(rgk1/P1*PQdif1,1,sum))+sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum))+sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum))
        if (NonUniform==T){
          Gamsco = c(sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum)))
        }
        if (m==2){
          Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        }else
          Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco,Gamsco,Betsco)# ,Betsco
          grad <- c(Dsco,Asco,Gamsco,Betsco)# ,Betsco
        } else {
          minusgrad <- -c(Dsco,Asco, Betsco)
          grad <- c(Dsco,Asco, Betsco)
        }
        
        if (NonUniform==T){
          FI <- matrix(0,length(minusgrad),length(minusgrad))#FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
          if (m==2){
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
              FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
            }
          } else {
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
              for (nn in m:(m+r-1)){
                FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
              }
            }
            if (m>2){
              for (mm in 2:(m-1)){
                FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
              }
              for (mm in 1:(m-2)){
                FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
              }
            }
          }
          FI[m,m] <- -sum(ng1*X[,1]*X[,1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
          FI[3,3]=FI[3,2]=FI[2,3]=-sum(ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum))
          FI[4,4]=FI[4,2]=FI[2,4]=-sum(ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
          FI[3,1]=FI[1,3]=sum(ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1]))
          FI[4,1]=FI[1,4]=sum(ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
          for (mm in 5:6){
            #(5,5) (6,6)
            FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
            #(2,5) (5,2) (2,6) (6,2)
            FI[m,mm] <- FI[mm,m] <- FI[mm,mm-2] <- FI[mm-2,mm] <- c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-4]
            #(1,5) (5,1) (1,6) (6,1)
            FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
          }
        } else {
          FI <- matrix(0,m+(m-1)*2,m+(m-1)*2)
          if (m==2){
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
              FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
            }
            for (mm in (m+1):(m+(m-1)*2)){
              FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
              FI[m,mm] <- FI[mm,m] <- c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-m]
              FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
            }
          } else {
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
              for (nn in m:(m+r-1)){
                FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
              }
            }
            if (m>2){
              for (mm in 2:(m-1)){
                FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
              }
              for (mm in 1:(m-2)){
                FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
              }
            }
          }
          FI[m,m] <- -sum(ng1*X[,1]*X[,1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
        }
        #sqrt(sum((grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6])^2))
        #gp=grad-c(d,a,gam[,1],bet)%*%FI
        gp=grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6]
        l2gp=sqrt(innerProduct(grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6],grad[3:6]-c(gam[,1],bet)%*%FI[3:6,3:6]))
        if (l2gp<eta){
          gam[,1]=bet[1:2]=0
          add <- qr.solve(FI,minusgrad)
          add[3:6]=0
          d=d+add[1:(m-1)]
          a=a+add[m]
        } else {
          add <- qr.solve(FI,minusgrad)
          direction=qr.solve(FI[3:6,3:6],as.vector(eta*gp/l2gp-grad[3:6]))
          #direction=zapsmall(c(c(gam[,1],bet),direction), digits = 9)[5:8]
          #for (ii in 1:4){
          #  if (abs(direction[ii])<1e-5){
          #    direction[ii]=0
          #  }
          #}
          linesch.alpha=1
          #if (direction[1] !=0 |direction[2] !=0 |direction[3] !=0 |direction[4] !=0 )
          #Armijo line search
          #{
          #  deltaj= direction%*%minusgrad[3:6]+eta*sqrt(innerProduct(c(gam[,1],bet)+direction,c(gam[,1],bet)+direction))-eta*sqrt(innerProduct(c(gam[,1],bet),c(gam[,1],bet)))
          #  steplength.vec=c(0.5^c(0:10))
          #  linesch=matrix(0,length(steplength.vec),4)
          #  for (ii in 1:length(steplength.vec)){
          #    linesch[ii,]=direction*steplength.vec[ii]
          #  }
          #  updtedtau.vec=matrix(rep(c(gam[,1],bet),length(steplength.vec)),length(steplength.vec),4,byrow=T)+linesch
          #  likelihood.vec=numeric(length(steplength.vec))
          #  for (ii in 1:length(steplength.vec)){
          #    sumoverk1=numeric(G)
          #    for(g in 1:G){
          #      sumoverk1[g]=rgk1[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1])))),0)))
          #    }
          #    sumoverk2=numeric(G)
          #    for(g in 1:G){
          #      sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(updtedtau.vec[ii,1]%*%X[g,1])+ updtedtau.vec[ii,3]))),0)))
          #    }
          #    sumoverk3=numeric(G)
          #    for(g in 1:G){
          #      sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(updtedtau.vec[ii,2]%*%X[g,1])+ updtedtau.vec[ii,4]))),0)))
          #    }
          #    likelihood.vec[ii]=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)
          #  }
          #  sumoverk2=numeric(G)
          #  for(g in 1:G){
          #    sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(gam[1,1]%*%X[g,1])+ bet[1]))),0)))
          #  }
          #  sumoverk3=numeric(G)
          #  for(g in 1:G){
          #    sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,1]) + rep(gam[2,1]%*%X[g,1])+ bet[2]))),0)))
          #  }
          #  lh0=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)
          #  Armijols=numeric(length(steplength.vec))
          #  for (ii in 1:length(steplength.vec)){
          #    Armijols[ii]=(lh0-likelihood.vec[ii]+eta*sqrt(innerProduct(c(gam[,1],bet)+direction*steplength.vec[ii],c(gam[,1],bet)+direction*steplength.vec[ii]))-eta*sqrt(innerProduct(c(gam[,1],bet),c(gam[,1],bet)))<=steplength.vec[ii]*0.1* deltaj)
          #  }
          #  if (max(Armijols)==0){
          #    linesch.alpha=0
          #  } else {
          #    linesch.alpha=max(steplength.vec[which(Armijols==1)])
          #  }
          #}
          d=d+add[1:(m-1)]
          a=a+add[m]
          add[3:6]=linesch.alpha*direction
          #add[3:6]=direction
          gam[,1]=gam[,1]+add[3:4]
          bet=bet+add[5:6]
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[1]
      gam=gam*Tau[1]
      
      gra[j,1] <- a
      grd[j,] <- d
      grgamma[,,j] <-gam
      grbeta[j,]=bet
    }
    
    # dimension2
    for (j in (r+J/2):J)
    {
      d <- grd[j,] 
      a <- gra[j,2]
      gam=grgamma[,,j]
      bet=grbeta[j,]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(k in 1:m){
        rLiA[,,k]=Xijk[,j,k]*LiA/Pi
      }
      rgk <- apply(rLiA,c(2,3),sum)
      rgk1 <- apply(rLiA[1:N1,,],c(2,3),sum)
      rgk2 <- apply(rLiA[(N1+1):(N1+N2),,],c(2,3),sum)
      rgk3 <- apply(rLiA[(N1+N2+1):(N1+N2+N3),,],c(2,3),sum)
      
      #M-step loop starts for item j
      miter <- 0
      add <- max.tol+1
      while(sum(abs(add))>max.tol)
      {
        miter <- miter+1
        for(g in 1:G){
          Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y1%*%bet+rep((y1%*%gam)%*%X[g,]))))
          P1[g,] <- -diff(c(1,Pstar1[g,],0))
          Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y2%*%bet+rep((y2%*%gam)%*%X[g,]))))
          P2[g,] <- -diff(c(1,Pstar2[g,],0))
          Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y3%*%bet+rep((y3%*%gam)%*%X[g,]))))
          P3[g,] <- -diff(c(1,Pstar3[g,],0))
        }
        Qstar <- 1-Pstar
        Qstar1 <- 1-Pstar1
        Qstar2 <- 1-Pstar2
        Qstar3 <- 1-Pstar3
        
        
        #calculating the score vector
        if (m==2){
          Dsco <- sum(apply(rgk1/P1,1,diff)*Pstar1*Qstar1+apply(rgk2/P2,1,diff)*Pstar2*Qstar2+apply(rgk3/P3,1,diff)*Pstar3*Qstar3)
        }else
          Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        PQdif1 <- -t(apply(cbind(0,Pstar1*Qstar1,0),1, diff))
        PQdif2 <- -t(apply(cbind(0,Pstar2*Qstar2,0),1, diff))
        PQdif3 <- -t(apply(cbind(0,Pstar3*Qstar3,0),1, diff))
        Asco = sum(X[,2]*apply(rgk1/P1*PQdif1,1,sum))+sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum))+sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum))
        if (NonUniform==T){
          Gamsco = c(sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum)))
        }
        if (m==2){
          Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        }else
          Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)#, Betsco
          grad <- c(Dsco,Asco,Gamsco,Betsco)#, Betsco
        } else {
          minusgrad <- -c(Dsco,Asco, Betsco)
          grad <- c(Dsco,Asco, Betsco)
        }
        
        if (NonUniform==T){
          FI <- matrix(0,length(minusgrad),length(minusgrad))#FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
          if (m==2){
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
              FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,2]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
            }
          } else { 
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
              for (nn in m:(m+r-1)){
                FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
              }
            }
            if (m>2){
              for (mm in 2:(m-1)){
                FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
              }
              for (mm in 1:(m-2)){
                FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
              }
            }
          }
          FI[m,m] <- -sum(ng1*X[,2]*X[,2]*apply(PQdif1^2/P1,1,sum)+ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum)+ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
          FI[3,3]=FI[3,2]=FI[2,3]=-sum(ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum))
          FI[4,4]=FI[4,2]=FI[2,4]=-sum(ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
          FI[3,1]=FI[1,3]=sum(ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1]))
          FI[4,1]=FI[1,4]=sum(ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
          for (mm in 5:6){
            #(5,5) (6,6)
            FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
            #(2,5) (5,2) (2,6) (6,2)
            FI[m,mm] <- FI[mm,m] <- FI[mm,mm-2] <- FI[mm-2,mm] <- c(sum(ng2*X[,2]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,2]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-4]
            #(1,5) (5,1) (1,6) (6,1)
            FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-4]
          }
        } else {
          FI <- matrix(0,m+(m-1)*2,m+(m-1)*2)
          if (m==2){
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
              FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,2]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
            }
            for (mm in (m+1):(m+(m-1)*2)){
              FI[mm,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
              FI[m,mm] <- FI[mm,m] <- c(sum(ng2*X[,2]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,2]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2])))[mm-m]
              FI[mm,m-1] <- FI[m-1,mm] <- c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2])))[mm-m]
            }
          } else {
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
              for (nn in m:(m+r-1)){
                FI[nn,mm] <- FI[mm,nn] <- sum(ng*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
              }
            }
            if (m>2){
              for (mm in 2:(m-1)){
                FI[mm,mm-1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
              }
              for (mm in 1:(m-2)){
                FI[mm,mm+1] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
              }
            }
          }
          FI[m,m] <- -sum(ng1*X[,2]*X[,2]*apply(PQdif1^2/P1,1,sum)+ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum)+ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
        }
        
        gp=grad[3:6]-c(gam[,2],bet)%*%FI[3:6,3:6]
        l2gp=sqrt(innerProduct(grad[3:6]-c(gam[,2],bet)%*%FI[3:6,3:6],grad[3:6]-c(gam[,2],bet)%*%FI[3:6,3:6]))
        if (l2gp<eta){
          gam[,2]=bet[1:2]=0
          add <- qr.solve(FI,minusgrad)
          add[3:6]=0
          d=d+add[1:(m-1)]
          a=a+add[m]
        } else {
          add <- qr.solve(FI,minusgrad)
          direction=qr.solve(FI[3:6,3:6],as.vector(eta*gp/l2gp-grad[3:6]))
          #direction=zapsmall(c(c(gam[,2],bet),direction), digits = 9)[5:8]
          #for (ii in 1:4){
          #  if (abs(direction[ii])<1e-5){
          #    direction[ii]=0
          #  }
          #}
          linesch.alpha=1
          #if (direction[1] !=0 |direction[2] !=0 |direction[3] !=0 |direction[4] !=0 )
          #{
          #  deltaj= direction%*%minusgrad[3:6]+eta*sqrt(innerProduct(c(gam[,2],bet)+direction,c(gam[,2],bet)+direction))-eta*sqrt(innerProduct(c(gam[,2],bet),c(gam[,2],bet)))
          #  steplength.vec=c(0.5^c(0:10))
          #  linesch=matrix(0,length(steplength.vec),4)
          #  for (ii in 1:length(steplength.vec)){
          #    linesch[ii,]=direction*steplength.vec[ii]
          #  }
          #  updtedtau.vec=matrix(rep(c(gam[,2],bet),length(steplength.vec)),length(steplength.vec),4,byrow=T)+linesch
          #  likelihood.vec=numeric(length(steplength.vec))
          #  for (ii in 1:length(steplength.vec)){
          #    sumoverk2=numeric(G)
          #    for(g in 1:G){
          #      sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(updtedtau.vec[ii,1]%*%X[g,2])+ updtedtau.vec[ii,3]))),0)))
          #    }
          #    sumoverk3=numeric(G)
          #    for(g in 1:G){
          #      sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(updtedtau.vec[ii,2]%*%X[g,2])+ updtedtau.vec[ii,4]))),0)))
          #    }
          #    likelihood.vec[ii]=sum(sumoverk2)+sum(sumoverk3)
          #  }
          #  sumoverk2=numeric(G)
          #  for(g in 1:G){
          #    sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(gam[1,2]%*%X[g,2])+ bet[1]))),0)))
          #  }
          #  sumoverk3=numeric(G)
          #  for(g in 1:G){
          #    sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,2]) + rep(gam[2,2]%*%X[g,2])+ bet[2]))),0)))
          #  }
          #  lh0=sum(sumoverk2)+sum(sumoverk3)
          #  Armijols=numeric(length(steplength.vec))
          #  for (ii in 1:length(steplength.vec)){
          #    Armijols[ii]=(lh0-likelihood.vec[ii]+eta*sqrt(innerProduct(c(gam[,2],bet)+direction*steplength.vec[ii],c(gam[,2],bet)+direction*steplength.vec[ii]))-eta*sqrt(innerProduct(c(gam[,2],bet),c(gam[,2],bet)))<=steplength.vec[ii]*0.1* deltaj)
          #  }
          #  if (max(Armijols)==0){
          #    linesch.alpha=0
          #  } else {
          #    linesch.alpha=max(steplength.vec[which(Armijols==1)])
          #  }
          #}
          d=d+add[1:(m-1)]
          a=a+add[m]
          add[3:6]=linesch.alpha*direction
          #add[3:6]=direction
          gam[,2]=gam[,2]+add[3:4]
          bet=bet+add[5:6]
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[2]
      gam=gam*Tau[2]
      
      gra[j,2] <- a
      grd[j,] <- d
      grgamma[,,j] <-gam
      grbeta[j,]=bet
    }
    
    
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.gamma <- abs(gammaold-grgamma)
    df.beta = abs(betaold-grbeta)
    iter <- iter+1
  }
  
  ##########################
  ###     Refit model    ###
  ##########################
  #get the sparsity structure
  sparsity=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  }
  anchor=which(sparsity[,1]==0)
  #Re-estimate by mirt
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  cmodel= multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchor]))
  estmirt = coef(cmodel,simplify=T)
  gramirt = estmirt$G1$items[,1:r]
  grdmirt = matrix(estmirt$G1$items[,(r+1):(r+m-1)],20,m-1)
  grgammamirt = array(0,dim=c(r,r,J))
  grgammamirt[1,,] = t(estmirt$G2$items[,1:r]-estmirt$G1$items[,1:r])
  grgammamirt[2,,] = t(estmirt$G3$items[,1:r]-estmirt$G1$items[,1:r])
  grbetamirt = matrix(0,J,2)
  grbetamirt[,1]=estmirt$G2$items[,(r+1):(r+m-1)]-estmirt$G1$items[,(r+1):(r+m-1)]
  grbetamirt[,2]=estmirt$G3$items[,(r+1):(r+m-1)]-estmirt$G1$items[,(r+1):(r+m-1)]
  Mu.gp1=estmirt$G1$means
  Mu.gp2=estmirt$G2$means
  Mu.gp3=estmirt$G3$means
  Sig.gp1=estmirt$G1$cov
  Sig.gp2=estmirt$G2$cov
  Sig.gp3=estmirt$G3$cov
  ######################
  ###       BIC      ###
  ######################
  est.d=grdmirt
  est.a=gramirt
  est.gamma=grgammamirt
  est.beta=grbetamirt
  
  if(min(resp)==0)
  {
    resp <- resp+1
  }
  G=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X <- matrix(gh[idx,1],nrow(idx),r)
  A1=dmvnorm(X,Mu.gp1,Sig.gp1)
  A2=dmvnorm(X,Mu.gp2,Sig.gp2)
  A3=dmvnorm(X,Mu.gp3,Sig.gp3)
  
  
  ##compute ng
  ng <-  numeric(G)
  #calculation of n_g 
  axmat=est.a%*%t(X) #a%*%X
  ygam1=ygam2=ygam3=matrix(0,20,2)
  for (j in 1:20){
    ygam1[j,]=y1%*%est.gamma[,,j]
  }
  for (j in 1:20){
    ygam2[j,]=y2%*%est.gamma[,,j]
  }
  for (j in 1:20){
    ygam3[j,]=y3%*%est.gamma[,,j]
  }
  ygamat1=ygam1%*%t(X)
  ygamat2=ygam2%*%t(X)
  ygamat3=ygam3%*%t(X)
  grbeta1=t(y1%*%t(est.beta))
  grbeta2=t(y2%*%t(est.beta))
  grbeta3=t(y3%*%t(est.beta))
  pstar1=pstar2=pstar3=array(double(J*(m-1)*G),dim = c(J,(m-1),G))
  p1=p2=p3=array(double(J*m*G),dim = c(J,m,G))
  for (g in 1:G)
  {
    pstar1[,,g] = 1/(1+exp(-(est.d+axmat[,g]+ygamat1[,g]+grbeta1)))
    p1[,,g] = t(-diff(rbind(rep(1,J),t(pstar1[,,g]),rep(0,J))))
  }
  for (g in 1:G)
  {
    pstar2[,,g] = 1/(1+exp(-(est.d+axmat[,g]+ygamat2[,g]+grbeta2)))
    p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
  }
  for (g in 1:G)
  {
    pstar3[,,g] = 1/(1+exp(-(est.d+axmat[,g]+ygamat3[,g]+grbeta3)))
    p3[,,g] = t(-diff(rbind(rep(1,J),t(pstar3[,,g]),rep(0,J))))
  }
  pij=array(double(J*G*N1),dim=c(J,G,N1))
  LiA=matrix(double(N*G),N,G)
  
  for (j in 1:J)
  {
    for (g in 1:G)
    {
      pij[j,g,]=p1[j,(resp[1:N1,j]),g]
    }
  }
  LiA[1:N1,]=t(apply(pij, c(2,3), prod)*A1)
  
  for (j in 1:J)
  {
    for (g in 1:G)
    {
      pij[j,g,]=p2[j,(resp[(N1+1):(N1+N2),j]),g]
    }
  }
  LiA[(N1+1):(N1+N2),]=t(apply(pij, c(2,3), prod)*A2)
  
  for (j in 1:J)
  {
    for (g in 1:G)
    {
      pij[j,g,]=p3[j,(resp[(N1+N2+1):(N1+N2+N3),j]),g]
    }
  }
  LiA[(N1+N2+1):(N1+N2+N3),]=t(apply(pij, c(2,3), prod)*A3)
  
  Pi = apply(LiA,1,sum)
  ng = apply(LiA/Pi,2,sum)
  ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
  ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
  ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
  
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for(j in 1:J)
  {
    d <- est.d[j,]
    a <- est.a[j,]
    gam=est.gamma[,,j]
    beta1=grbeta1[j,]
    beta2=grbeta2[j,]
    beta3=grbeta3[j,]
    
    rLiA <- array(double(N*G*m),dim = c(N,G,m))
    for(k in 1:m){
      rLiA[,,k]=Xijk[,j,k]*LiA/Pi
    }
    rgk <- apply(rLiA,c(2,3),sum)
    rgk1 <- apply(rLiA[1:N1,,],c(2,3),sum)
    rgk2 <- apply(rLiA[(N1+1):(N1+N2),,],c(2,3),sum)
    rgk3 <- apply(rLiA[(N1+N2+1):(N1+N2+N3),,],c(2,3),sum)
    
    sumoverk1=numeric(G)
    for(g in 1:G){
      sumoverk1[g]=rgk1[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,])+ beta1))),0)))
    }
    sumoverk2=numeric(G)
    for(g in 1:G){
      sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,])+rep(gam[1,]%*%X[g,])+ beta2))),0)))
    }
    sumoverk3=numeric(G)
    for(g in 1:G){
      sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,])+rep(gam[2,]%*%X[g,])+ beta3))),0)))
    }
    temp=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm=numeric(J) 
  for (j in 1:20){
    l0norm[j]=ifelse(apply(sparsity,1,sum)[j]==0,0,1)
  }
  dfl=0
  for (j in 3:J){
    dfl=dfl+3*(sqrt(sum((grgamma^2)[,,j])+sum((grbeta^2)[j,]))/sqrt(sum((grgamma000^2)[,,j])+sum((grbeta000^2)[j,])))
  }
  df=J+N+(sum(l0norm)+dfl)-1
  BIC=-2*sum(lh)+df*log(N*J)
  Bias=c(colSums(est.a-Amat1)/10,colMeans(est.d-Dmat1))
  RMSE=c(sqrt(colSums((est.a-Amat1)^2)/10),sqrt(colMeans((est.d-Dmat1)^2)))
  
  #output esimates and number of iterations
  return(list(est=cbind(est.a,est.d),Gamma=est.gamma,Beta=est.beta,mean1=Mu.gp1,mean2=Mu.gp2,mean3=Mu.gp3,Corr1=Sig.gp1,Corr2=Sig.gp2,Corr3=Sig.gp3,iter=iter,bic=BIC,bias=Bias,RMSE=RMSE))
}
#end of function

s <- 'F1 = 1,3-11
          F2 = 2,12-20
          COV = F1*F2'
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
md00 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:r]))
gra000=coef(md00,simplify=T)$G1$items[,c("a1","a2")]
grd000=matrix(coef(md00,simplify=T)$G1$items[,c("d")],20,1)
grgamma000=array(0,dim=c(r,r,J))
grgamma000[1,1,3:11]=(coef(md00,simplify=T)$G2$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[3:11,1]
grgamma000[2,1,3:11]=(coef(md00,simplify=T)$G3$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[3:11,1]
grgamma000[1,2,12:20]=(coef(md00,simplify=T)$G2$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[12:20,2]
grgamma000[2,2,12:20]=(coef(md00,simplify=T)$G3$items[,c("a1","a2")]-coef(md00,simplify=T)$G1$items[,c("a1","a2")])[12:20,2]
grbeta000=matrix(0,J,2)
grbeta000[,1]=matrix(coef(md00,simplify=T)$G2$items[,c("d")],20,1)-matrix(coef(md00,simplify=T)$G1$items[,c("d")],20,1)
grbeta000[,2]=matrix(coef(md00,simplify=T)$G3$items[,c("d")],20,1)-matrix(coef(md00,simplify=T)$G1$items[,c("d")],20,1)
Sig100=coef(md00,simplify=T)$G1$cov
Sig200=coef(md00,simplify=T)$G2$cov
Sig300=coef(md00,simplify=T)$G3$cov
mu100=coef(md00,simplify=T)$G1$means
mu200=coef(md00,simplify=T)$G2$means
mu300=coef(md00,simplify=T)$G3$means

for (rep in 2:50){
  set.seed(rep*100)
  Theta=mvrnorm(n=N,mu=c(0,0),Sigma = matrix(c(1,0.85,0.85,1),2,2))
  Theta1=Theta[1:N1,]
  Theta2=Theta[(N1+1):(N1+N2),]
  Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
  datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
  datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
  datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
  resp=rbind(datasetR,datasetF1,datasetF2)
  
  r=2
  m=2
  eta.vec=seq(21,48,3)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*m*length(eta.vec)),dim = c(J,m,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    sim=ipest1(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=T,gra00=gra00,grd00=grd00,grgamma00=grgamma00,grbeta00=grbeta00,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
    bics[k]=sim$bic
    ADmat[,,k]=sim$est
    Gammas[,,,k]=sim$Gamma
    Betas[,,k]=sim$Beta
    biass[k,]=sim$bias
    RMSEs[k,]=sim$RMSE
  }
  
  
  
  kk=which.min(bics)
  print(eta.vec[kk])
  print(Betas[,,kk])
  
  write.csv(eta.vec[kk],file = paste("eta9_",rep))
  write.csv(ADmat[,,kk],file = paste("ADmat9_",rep))
  write.csv(Gammas[,,,kk],file = paste("Gamma9_",rep))
  write.csv(Betas[,,kk],file = paste("Beta9_",rep))
  
  eta.13[rep]=eta.vec[kk]
  ADmat.13[,,rep]=ADmat[,,kk]
  Gammas.13[,,,rep]=Gammas[,,,kk]
  Betas.13[,,rep]=Betas[,,kk]
  biass.13[rep,]=biass[kk,]
  RMSEs.13[rep,]=RMSEs[kk,]
}
