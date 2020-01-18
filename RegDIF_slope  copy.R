# Dataset #3 (2 Non-uniform DIF items per scale)
set.seed(1)
A11=c(runif(1,1.5,2.5),0,runif(9,1.5,2.5),rep(0,9))
A12=A13=A11
A21=c(0,runif(1,1.5,2.5),rep(0,9),runif(9,1.5,2.5))
A22=A23=A21
D11=rnorm(J,0,1)
D12=D13=D11
# focal group has smaller slopes and more extreme intercepts
A12[5]=A11[5]- 0.3
A22[12]=A21[12]- 0.3
A22[13]=A21[13]- 0.3
A13[4]=A11[4]- 0.6
A13[5]=A11[5]- 0.6
A23[12]=A21[12]- 0.6
A23[13]=A21[13]- 0.6
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

##################################################
### Lasso DIF Quadrature (mirt start, re-est) ####
##################################################

library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(irtoys)

soft=function(s, tau) {      
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

#starting values by mirt
s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
cmodel <- mirt.model(s)
md0 <- multipleGroup(resp, cmodel, group = Group,invariance=c('slopes',colnames(resp)[1:r]))
gra0=coef(md0,simplify=T)$R$items[,1:r]
grd0=matrix(coef(md0,simplify=T)$R$items[,3],J,1)
grbeta0[,1]=matrix(coef(md0,simplify=T)$F1$items[,3]-coef(md0,simplify=T)$R$items[,3],J,1)
grbeta0[,2]=matrix(coef(md0,simplify=T)$F2$items[,3]-coef(md0,simplify=T)$R$items[,3],J,1)
Sig0=coef(md0,simplify=T)$R$cov

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
grgamma0[1,1,3:11]=-0.5
grgamma0[2,1,3:11]=-1
grgamma0[1,2,12:20]=-0.5
grgamma0[2,2,12:20]=-1
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
ng <-  numeric(G)
#gra = azero
#grd = matrix(0,J,1)
gra=gra0
grd=grd0
#grbeta=grbeta0
Sig=Sig0
grgamma=grgamma0
grbeta=matrix(0,J,2)
#starting values for the MC samples form the prior
Mu=rep(0,r)
#X=theta #for check later

Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
df.a <- df.d  <- df.gamma <- df.Sig <- 1

iter <- 0

while(max(df.a)>eps & max(df.d)>eps&  max(df.gamma)>eps & max(df.Sig)>eps)
{
  aold <- gra
  dold <- grd
  #gammaold=grgamma
  gammaold=grgamma
  Sigold <- Sig
  A=dmvnorm(X,rep(0,r),Sig)
  
  #calculation of n_g 
  axmat=gra%*%t(X) #a%*%X
  ygam1=ygam2=ygam3=matrix(0,20,2)
  for (j in 1:20){
    ygam1[j,]=y1%*%grgamma[,,j]
  }
  for (j in 1:20){
    ygam2[j,]=y2%*%grgamma[,,j]
  }
  for (j in 1:20){
    ygam3[j,]=y3%*%grgamma[,,j]
  }
  ygamat1=ygam1%*%t(X)
  ygamat2=ygam2%*%t(X)
  ygamat3=ygam3%*%t(X)
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
  pij=matrix(double(J*G),J,G)
  LiA=matrix(double(N*G),N,G)
  
  for (i in 1:N1)
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p1[j,respi[j],g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A
  }
  
  for (i in (N1+1):(N1+N2))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p2[j,respi[j],g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A
  }
  
  for (i in (N1+N2+1):(N1+N2+N3))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p3[j,respi[j],g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A
  }
  
  Pi = apply(LiA,1,sum)
  ng = apply(LiA/Pi,2,sum)
  ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
  ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
  ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
  
  #update mu hat and Sigma hat
  Sigg=matrix(0,r,r)
  for (g in 1:G){
    Sigg=Sigg+X[g,]%*%t(X[g,])*ng[g]
  }
  Sig.hat=Sigg/N
  
  #scale 
  #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
  Tau=sqrt(diag(Sig.hat))
  Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
  #q_g_star
  #X=(X-mu.hat.mat)/Tau.mat
  Xstar=X/Tau.mat
  
  Sigg2=matrix(0,r,r)
  for (g in 1:G){
    Sigg2=Sigg2+Xstar[g,]%*%t(Xstar[g,])*ng[g]
  }
  Sig=Sigg2/N
  
  ##Constraint part
  #calculation of r_jgk
  for (j in 1:r)
  {
    d <- grd[j,]
    a <- gra[j,j]
    rLiA <- array(double(N*G*m),dim = c(N,G,m))
    for(i in 1:N){
      for (g in 1:G){
        rLiA[i,g,resp[i,j]] <- LiA[i,g]
      }
      rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
    for(i in 1:N){
      for (g in 1:G){
        rLiA[i,g,resp[i,j]] <- LiA[i,g]
      }
      rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
      #if (m==2){
      #  Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      #}else
      #  Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      if (NonUniform==T){
        minusgrad <- -c(Dsco,Asco,Gamsco)# ,Betsco
        grad <- c(Dsco,Asco,Gamsco)# ,Betsco
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
      
      
      add <- qr.solve(FI,minusgrad)
      d=d+add[1:(m-1)]
      a=a+add[m]
      gam0=gam[,1]+add[(m+r-1):(m+r)]
      for (mm in (m+r-1):(m+r)){
        add[mm]=soft(gam0[mm-m],-eta/FI[mm,mm])-gam[mm-m,1]
      }
      for (mm in (m+r-1):(m+r)){
        gam[mm-m,1]=soft(gam0[mm-m],-eta/FI[mm,mm])
      }
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[1]
    gam=gam*Tau[1]
    
    gra[j,1] <- a
    grd[j,] <- d
    grgamma[,,j] <-gam
  }
  
  # dimension2
  for (j in (r+J/2):J)
  {
    d <- grd[j,] 
    a <- gra[j,2]
    gam=grgamma[,,j]
    bet=grbeta[j,]
    rLiA <- array(double(N*G*m),dim = c(N,G,m))
    for(i in 1:N){
      for (g in 1:G){
        rLiA[i,g,resp[i,j]] <- LiA[i,g]
      }
      rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
      #if (m==2){
      #  Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      #}else
      #  Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      if (NonUniform==T){
        minusgrad <- -c(Dsco,Asco, Gamsco)#, Betsco
        grad <- c(Dsco,Asco,Gamsco)#, Betsco
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
      
      
      add <- qr.solve(FI,minusgrad)
      d=d+add[1:(m-1)]
      a=a+add[m]
      gam0=gam[,2]+add[(m+r-1):(m+r)]
      for (mm in (m+r-1):(m+r)){
        add[mm]=soft(gam0[mm-m],-eta/FI[mm,mm])-gam[mm-m,2]
      }
      for (mm in (m+r-1):(m+r)){
        gam[mm-m,2]=soft(gam0[mm-m],-eta/FI[mm,mm])
      }
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[2]
    gam=gam*Tau[2]
    
    gra[j,2] <- a
    grd[j,] <- d
    grgamma[,,j] <-gam
  }
  
  
  df.d <- abs(dold-grd)
  df.a <- abs(aold-gra)
  df.gamma <- abs(gammaold-grgamma)
  df.Sig = abs(Sigold-Sig)
  iter <- iter+1
}


gra00=gra
grd00=grd
grgamma00=grgamma
Sig00=Sig

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

ipest1 <- function(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=gra00,grd00=grd00,grgamma00=grgamma00,Sig00=Sig00)
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
  ng <-  numeric(G)
  #gra = azero
  #grd = matrix(0,J,1)
  gra=gra00
  grd=grd00
  grgamma=grgamma00
  grbeta=matrix(0,J,2)
  Sig=Sig00
  
  #starting values for the MC samples form the prior
  Mu=rep(0,r)
  #X=theta #for check later
  
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps & max(df.d)>eps & max(df.gamma)>eps & max(df.Sig)>eps)
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    Sigold <- Sig
    A=dmvnorm(X,rep(0,r),Sig)
    
    #calculation of n_g 
    axmat=gra%*%t(X) #a%*%X
    ygam1=ygam2=ygam3=matrix(0,20,2)
    for (j in 1:20){
      ygam1[j,]=y1%*%grgamma[,,j]
    }
    for (j in 1:20){
      ygam2[j,]=y2%*%grgamma[,,j]
    }
    for (j in 1:20){
      ygam3[j,]=y3%*%grgamma[,,j]
    }
    ygamat1=ygam1%*%t(X)
    ygamat2=ygam2%*%t(X)
    ygamat3=ygam3%*%t(X)
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
    pij=matrix(double(J*G),J,G)
    LiA=matrix(double(N*G),N,G)
    
    for (i in 1:N1)
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p1[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    
    for (i in (N1+1):(N1+N2))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p2[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    
    for (i in (N1+N2+1):(N1+N2+N3))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p3[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    
    Pi = apply(LiA,1,sum)
    ng = apply(LiA/Pi,2,sum)
    ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
    ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
    ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
    
    #update mu hat and Sigma hat
    Sigg=matrix(0,r,r)
    for (g in 1:G){
      Sigg=Sigg+X[g,]%*%t(X[g,])*ng[g]
    }
    Sig.hat=Sigg/N
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    
    Sigg2=matrix(0,r,r)
    for (g in 1:G){
      Sigg2=Sigg2+Xstar[g,]%*%t(Xstar[g,])*ng[g]
    }
    Sig=Sigg2/N
    
    ##Constraint part
    #calculation of r_jgk
    for (j in 1:r)
    {
      d <- grd[j,]
      a <- gra[j,j]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
        #if (m==2){
        #  Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        #}else
        #  Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco,Gamsco)# ,Betsco
          grad <- c(Dsco,Asco,Gamsco)# ,Betsco
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
        
        
        add <- qr.solve(FI,minusgrad)
        d=d+add[1:(m-1)]
        a=a+add[m]
        gam0=gam[,1]+add[(m+r-1):(m+r)]
        for (mm in (m+r-1):(m+r)){
          add[mm]=soft(gam0[mm-m],-eta/FI[mm,mm])-gam[mm-m,1]
        }
        for (mm in (m+r-1):(m+r)){
          gam[mm-m,1]=soft(gam0[mm-m],-eta/FI[mm,mm])
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[1]
      gam=gam*Tau[1]
      
      gra[j,1] <- a
      grd[j,] <- d
      grgamma[,,j] <-gam
    }
    
    # dimension2
    for (j in (r+J/2):J)
    {
      d <- grd[j,] 
      a <- gra[j,2]
      gam=grgamma[,,j]
      bet=grbeta[j,]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
        #if (m==2){
        #  Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        #}else
        #  Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco, Gamsco)#, Betsco
          grad <- c(Dsco,Asco,Gamsco)#, Betsco
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
        
        
        add <- qr.solve(FI,minusgrad)
        d=d+add[1:(m-1)]
        a=a+add[m]
        gam0=gam[,2]+add[(m+r-1):(m+r)]
        for (mm in (m+r-1):(m+r)){
          add[mm]=soft(gam0[mm-m],-eta/FI[mm,mm])-gam[mm-m,2]
        }
        for (mm in (m+r-1):(m+r)){
          gam[mm-m,2]=soft(gam0[mm-m],-eta/FI[mm,mm])
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[2]
      gam=gam*Tau[2]
      
      gra[j,2] <- a
      grd[j,] <- d
      grgamma[,,j] <-gam
    }
    
    
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.gamma <- abs(gammaold-grgamma)
    df.Sig = abs(Sigold-Sig)
    iter <- iter+1
  }
  
  ###############
  # Refit model #
  ###############
  
  #get the sparsity structure
  sparsity=grgamma
  for (j in 1:J){
    for (rr in 1:2){
      for (nn in 1:2){
        sparsity[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.gamma  <- df.Sig <- 1
  
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps & max(df.d)>eps & max(df.gamma)>eps & max(df.Sig)>eps)
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    Sigold <- Sig
    A=dmvnorm(X,rep(0,r),Sig)
    
    #calculation of n_g 
    axmat=gra%*%t(X) #a%*%X
    ygam1=ygam2=ygam3=matrix(0,20,2)
    for (j in 1:20){
      ygam1[j,]=y1%*%grgamma[,,j]
    }
    for (j in 1:20){
      ygam2[j,]=y2%*%grgamma[,,j]
    }
    for (j in 1:20){
      ygam3[j,]=y3%*%grgamma[,,j]
    }
    ygamat1=ygam1%*%t(X)
    ygamat2=ygam2%*%t(X)
    ygamat3=ygam3%*%t(X)
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
    pij=matrix(double(J*G),J,G)
    LiA=matrix(double(N*G),N,G)
    
    for (i in 1:N1)
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p1[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    
    for (i in (N1+1):(N1+N2))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p2[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    
    for (i in (N1+N2+1):(N1+N2+N3))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p3[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    
    Pi = apply(LiA,1,sum)
    ng = apply(LiA/Pi,2,sum)
    ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
    ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
    ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
    
    #update mu hat and Sigma hat
    Sigg=matrix(0,r,r)
    for (g in 1:G){
      Sigg=Sigg+X[g,]%*%t(X[g,])*ng[g]
    }
    Sig.hat=Sigg/N
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    
    Sigg2=matrix(0,r,r)
    for (g in 1:G){
      Sigg2=Sigg2+Xstar[g,]%*%t(Xstar[g,])*ng[g]
    }
    Sig=Sigg2/N
    
    ##Constraint part
    #calculation of r_jgk
    for (j in 1:r)
    {
      d <- grd[j,]
      a <- gra[j,j]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
      l0normgamj=0
      for(l in 1:2) 
      {
        for (n in 1:2){
          l0normgamj=l0normgamj+(sparsity[n,l,j]!=0)
        }
      }
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
        Qstar1 <- 1-Pstar1
        Qstar2 <- 1-Pstar2
        Qstar3 <- 1-Pstar3
        
        
        #calculating the score vector
        if (m==2){
          Dsco <- sum(apply(rgk1/P1,1,diff)*Pstar1*Qstar1+apply(rgk2/P2,1,diff)*Pstar2*Qstar2+apply(rgk3/P3,1,diff)*Pstar3*Qstar3)
        }else
          Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        PQdif1 <- -t(apply(cbind(0,Pstar1*Qstar1,0),1, diff))
        PQdif2 <- -t(apply(cbind(0,Pstar2*Qstar2,0),1, diff))
        PQdif3 <- -t(apply(cbind(0,Pstar3*Qstar3,0),1, diff))
        Asco = sum(X[,1]*apply(rgk1/P1*PQdif1,1,sum))+sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum))+sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum))
        if (NonUniform==T){
          if (l0normgamj==1){
            Gamsco = sparsity[,1,j]%*%c(sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum)))
          }else if (l0normgamj==2){
            Gamsco = c(sum(X[,1]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,1]*apply(rgk3/P3*PQdif3,1,sum)))
          } else{
            Gamsco =numeric(0)
          }
        }
        #if (m==2){
        #  Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        #}else
        #  Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco,Gamsco)# ,Betsco
          grad <- c(Dsco,Asco,Gamsco)# ,Betsco
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
          if (l0normgamj==1){
            FI[3,3]=FI[3,2]=FI[2,3]=sparsity[,1,j]%*%c(-sum(ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum)),-sum(ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum)))
            FI[3,1]=FI[1,3]=sparsity[,1,j]%*%c(sum(ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])),sum(ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1])))
          }
          if (l0normgamj==2){
            FI[3,3]=FI[3,2]=FI[2,3]=-sum(ng2*X[,1]*X[,1]*apply(PQdif2^2/P2,1,sum))
            FI[4,4]=FI[4,2]=FI[2,4]=-sum(ng3*X[,1]*X[,1]*apply(PQdif3^2/P3,1,sum))
            FI[3,1]=FI[1,3]=sum(ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1]))
            FI[4,1]=FI[1,4]=sum(ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
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
        
        
        add <- qr.solve(FI,minusgrad)
        d=d+add[1:(m-1)]
        a=a+add[m]
        if (l0normgamj==2){
          gam[,1]=gam[,1]+add[3:4]
        } else if (l0normgamj==1){
          gam[which(sparsity[,1,j]==1),1]=gam[which(sparsity[,1,j]==1),1]+add[3]
        }
        
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[1]
      if (l0normgamj!=0){
        gam=gam*Tau[1]
      }
      gra[j,1] <- a
      grd[j,] <- d
      grgamma[,,j] <-gam
    }
    
    # dimension2
    for (j in (r+J/2):J)
    {
      d <- grd[j,] 
      a <- gra[j,2]
      gam=grgamma[,,j]
      bet=grbeta[j,]
      l0normgamj=0
      for(l in 1:2) 
      {
        for (n in 1:2){
          l0normgamj=l0normgamj+(sparsity[n,l,j]!=0)
        }
      }
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
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
        Qstar1 <- 1-Pstar1
        Qstar2 <- 1-Pstar2
        Qstar3 <- 1-Pstar3
        
        
        #calculating the score vector
        if (m==2){
          Dsco <- sum(apply(rgk1/P1,1,diff)*Pstar1*Qstar1+apply(rgk2/P2,1,diff)*Pstar2*Qstar2+apply(rgk3/P3,1,diff)*Pstar3*Qstar3)
        }else
          Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        PQdif1 <- -t(apply(cbind(0,Pstar1*Qstar1,0),1, diff))
        PQdif2 <- -t(apply(cbind(0,Pstar2*Qstar2,0),1, diff))
        PQdif3 <- -t(apply(cbind(0,Pstar3*Qstar3,0),1, diff))
        Asco = sum(X[,2]*apply(rgk1/P1*PQdif1,1,sum))+sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum))+sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum))
        if (NonUniform==T) {
          if (l0normgamj==1) {
            Gamsco = sparsity[,2,j]%*%c(sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum)))
          } else if (l0normgamj==2) {
            Gamsco = c(sum(X[,2]*apply(rgk2/P2*PQdif2,1,sum)),sum(X[,2]*apply(rgk3/P3*PQdif3,1,sum)))
          } else {
            Gamsco =numeric(0)
          }
        }
        #if (m==2){
        #  Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        #}else
        #  Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco, Gamsco)#, Betsco
          grad <- c(Dsco,Asco,Gamsco)#, Betsco
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
          if (l0normgamj==1){
            FI[3,3]=FI[3,2]=FI[2,3]=sparsity[,2,j]%*%c(-sum(ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum)),-sum(ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum)))
            FI[3,1]=FI[1,3]=sparsity[,2,j]%*%c(sum(ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])),sum(ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1])))
          }
          if (l0normgamj==2){
            FI[3,3]=FI[3,2]=FI[2,3]=-sum(ng2*X[,2]*X[,2]*apply(PQdif2^2/P2,1,sum))
            FI[4,4]=FI[4,2]=FI[2,4]=-sum(ng3*X[,2]*X[,2]*apply(PQdif3^2/P3,1,sum))
            FI[3,1]=FI[1,3]=sum(ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1]))
            FI[4,1]=FI[1,4]=sum(ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
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
        
        
        add <- qr.solve(FI,minusgrad)
        d=d+add[1:(m-1)]
        a=a+add[m]
        if (l0normgamj==2){
          gam[,2]=gam[,2]+add[3:4]
        } 
        if (l0normgamj==1){
          gam[which(sparsity[,2,j]==1),2]=gam[which(sparsity[,2,j]==1),2]+add[3]
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[2]
      if (l0normgamj!=0){
        gam=gam*Tau[2]
      }   
      gra[j,2] <- a
      grd[j,] <- d
      grgamma[,,j] <-gam
    }
    
    
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.gamma <- abs(gammaold-grgamma)
    df.Sig = abs(Sigold-Sig)
    iter <- iter+1
  }
  
  ######################
  ###       BIC      ###
  ######################
  est.d=grd
  est.a=gra
  est.gamma=grgamma
  est.beta=grbeta
  
  if(min(resp)==0)
  {
    resp <- resp+1
  }
  G=length(X1)^r
  gh=t(matrix(rep(X1,r),r,length(X1),byrow = T))
  idx <- as.matrix(expand.grid(rep(list(1:length(X1)),r)))
  X <- matrix(gh[idx,1],nrow(idx),r)
  A=dmvnorm(X,rep(0,r),Sig)
  
  
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
  pij=matrix(double(J*G),J,G)
  LiA=matrix(double(N*G),N,G)
  
  for (i in 1:N1)
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p1[j,respi[j],g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A
  }
  
  for (i in (N1+1):(N1+N2))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p2[j,respi[j],g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A
  }
  
  for (i in (N1+N2+1):(N1+N2+N3))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p3[j,respi[j],g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A
  }
  
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
    for(i in 1:N){
      for (g in 1:G){
        rLiA[i,g,resp[i,j]] <- LiA[i,g]
      }
      rLiA[i,,] <- rLiA[i,,]/Pi[i]
    }
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
  l0norm=0
  for(i in 1:J) 
  {
    for(j in 1:2)
    {
      for(k in 1:2){
        l0norm=l0norm+(est.gamma[j,k,i]!=0)
      }
    }
  }
  
  BIC=-2*sum(lh)+l0norm*log(N)
  
  Bias=c(colSums(gra-Amat1)/10,colMeans(grd-Dmat1))
  RMSE=c(sqrt(colSums((gra-Amat1)^2)/10),sqrt(colMeans((grd-Dmat1)^2)))
  
  #output esimates and number of iterations
  return(list(est=cbind(gra,grd),Gamma=grgamma,Corr=Sig,iter=iter,bic=BIC,bias=Bias,RMSE=RMSE))
}
#end of function

r=2
m=2
eta.vec=seq(15,29,2)
bics=rep(0,length(eta.vec))
Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
biass=matrix(0,length(eta.vec),3)
RMSEs=matrix(0,length(eta.vec),3)
for (k in 1:length(eta.vec))
{
  eta=eta.vec[k]
  sim=ipest1(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=T,gra00=gra00,grd00=grd00,grgamma00=grgamma00,Sig00=Sig00)
  bics[k]=sim$bic
  Gammas[,,,k]=sim$Gamma
  biass[k,]=sim$bias
  RMSEs[k,]=sim$RMSE
}

kk=which.min(bics)
eta.vec[kk]
Gammas[,,,kk]
biass[kk,]
RMSEs[kk,]






#########################################
#####                               #####
#####  Comparing with mirt LRT add  #####
#####                               #####
#########################################
s <- 'F1 = 1,3-11
          F2 = 2,12-20
          COV = F1*F2'
cmodel <- mirt.model(s)
mirtCluster()
md.noncons0 <- multipleGroup(resp, cmodel, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
DIF(md.noncons0, which.par = c('a1'), Wald=TRUE, p.adjust = 'fdr',scheme = 'add',items2test=c(3:11)) 
DIF(md.noncons0, which.par = c('a2'), Wald=TRUE, p.adjust = 'fdr',scheme = 'add',items2test=c(12:20)) 
mirt100[,8]


mirt.p.mat=matrix(0,16,18)
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
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  cmodel <- mirt.model(s)
  md.noncons0 <- multipleGroup(resp, cmodel, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:r]))
  mirt.p.mat[(i),]=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,8]
  
}
sum(mirt.p.mat[,c(2,3,10,11)]<0.01)/(16*4)
sum(mirt.p.mat[,-c(2,3,10,11)]<0.05)/(16*14)



mirt.p.mat=matrix(0,20,18) # 20 is the number of replications
bias.mirt=matrix(0,20,3)
rmse.mirt=matrix(0,20,3)
difrec.mirt.fn=matrix(0,20,3) # DIF magnitude recovery with false negative included
difrec.mirt=matrix(0,20,3) # DIF magnitude recovery without false negative (the results are same as false negative included case)
difmag.mirt.fp=matrix(0,20,2) # DIF magnitude recovery among those none DIF items (false positive)
for (i in 2:20){
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
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #md.noncons0 <- multipleGroup(resp, cmodel, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('intercepts',colnames(resp)[1:r]))
  mirt.p.mat[(i),1:9]=DIF(md.noncons0, which.par = c('a1'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:11))[,8]
  mirt.p.mat[(i),10:18]=DIF(md.noncons0, which.par = c('a2'), p.adjust = 'fdr',scheme = 'add',items2test=c(12:20))[,8]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('intercepts',colnames(resp)[-(which(mirt.p.mat[i,]<0.05)+2)]))
  bias.mirt[i,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt[i,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt[i,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt[i,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn[i,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[c(4,5)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[c(12,13)]-0.5)))
  difrec.mirt.fn[i,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[c(4,5)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[c(12,13)]-1)))
  difrec.mirt.fn[i,1]=mean(c(difrec.mirt.fn[i,2],difrec.mirt.fn[i,3]))
  #difmag.mirt.fp[i,1]= mean(abs((coef(md.refit0,simplify=T)$G1$items[,1:2]-coef(md.refit0,simplify=T)$G2$items[,1:2])[-c(4,5,12,13),]))
  #difmag.mirt.fp[i,2]= mean(abs((coef(md.refit0,simplify=T)$G1$items[,1:2]-coef(md.refit0,simplify=T)$G3$items[,1:2])[-c(4,5,12,13),]))
}
sum(mirt.p.mat[,c(2,3,10,11)]<0.05)/(20*4)
sum(mirt.p.mat[,-c(2,3,10,11)]<0.05)/(20*14)
colMeans(bias.mirt)
colMeans(rmse.mirt)
colMeans(difrec.mirt.fn)
colMeans(na.omit(difmag.mirt.fp))

