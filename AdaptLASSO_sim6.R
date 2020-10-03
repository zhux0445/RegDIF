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
responses=read.csv("RESP6.csv",row.names = 1)

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

ipest1 <- function(resp,m,r,eta,lam,eps =1e-3,max.tol=1e-7,NonUniform=T,gra00=gra00,grd00=grd00,grgamma00=grgamma00,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
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
  grbeta=matrix(0,J,2)
  Sig.gp1=Sig100; Sig.gp2=Sig200; Sig.gp3=Sig300
  
  #starting values for the MC samples form the prior
  Mu.gp1=mu100; Mu.gp2=mu200; Mu.gp3=mu300
  #X=theta #for check later
  
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.gamma <- df.beta <- 1
  
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps & max(df.d)>eps & max(df.gamma)>eps)
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    A1=dmvnorm(X,Mu.gp1,Sig.gp1)
    A2=dmvnorm(X,Mu.gp2,Sig.gp2)
    A3=dmvnorm(X,Mu.gp3,Sig.gp3)
    
    #calculation of n_g 
    axmat=eigenMapMatMult(gra,t(X)) #a%*%X
    ygam1=ygam2=ygam3=matrix(0,J,r)
    ygam2=t(grgamma[1,,])
    ygam3=t(grgamma[2,,])
    ygamat1=eigenMapMatMult(ygam1,t(X))
    ygamat2=eigenMapMatMult(ygam2,t(X))
    ygamat3=eigenMapMatMult(ygam3,t(X))
    grbeta1=matrix(0,J,m-1)
    grbeta2=matrix(0,J,m-1)
    grbeta3=matrix(0,J,m-1)
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
    
    #update Sigma hat
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
      gammle=grgamma00[,,j]
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
          add[mm]=soft(gam0[mm-m],-(eta/(abs(gammle[mm-m,1])^(lam)))/FI[mm,mm])-gam[mm-m,1]
        }
        for (mm in (m+r-1):(m+r)){
          gam[mm-m,1]=soft(gam0[mm-m],-(eta/(abs(gammle[mm-m,1])^(lam)))/FI[mm,mm])
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
      gammle=grgamma00[,,j]
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
          add[mm]=soft(gam0[mm-m],-(eta/(abs(gammle[mm-m,2])^(lam)))/FI[mm,mm])-gam[mm-m,2]
        }
        for (mm in (m+r-1):(m+r)){
          gam[mm-m,2]=soft(gam0[mm-m],-(eta/(abs(gammle[mm-m,2])^(lam)))/FI[mm,mm])
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
  df.a <- df.d  <- df.gamma  <- 1
  
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps & max(df.d)>eps & max(df.gamma)>eps)
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    A1=dmvnorm(X,Mu.gp1,Sig.gp1)
    A2=dmvnorm(X,Mu.gp2,Sig.gp2)
    A3=dmvnorm(X,Mu.gp3,Sig.gp3)
    
    #calculation of n_g 
    axmat=eigenMapMatMult(gra,t(X)) #a%*%X
    ygam1=ygam2=ygam3=matrix(0,J,r)
    ygam2=t(grgamma[1,,])
    ygam3=t(grgamma[2,,])
    ygamat1=eigenMapMatMult(ygam1,t(X))
    ygamat2=eigenMapMatMult(ygam2,t(X))
    ygamat3=eigenMapMatMult(ygam3,t(X))
    grbeta1=matrix(0,J,m-1)
    grbeta2=matrix(0,J,m-1)
    grbeta3=matrix(0,J,m-1)
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
    
    #update Sigma hat
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
      l0normgamj=0
      for(l in 1:2) 
      {
        for (n in 1:2){
          l0normgamj=l0normgamj+(sparsity[n,l,j]!=0)
        }
      }
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
    iter <- iter+1
  }
  
  ######################
  ###       GIC      ###
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
  A1=dmvnorm(X,Mu.gp1,Sig.gp1)
  A2=dmvnorm(X,Mu.gp2,Sig.gp2)
  A3=dmvnorm(X,Mu.gp3,Sig.gp3)
  
  
  ##compute ng
  ng <-  numeric(G)
  #calculation of n_g 
  axmat=eigenMapMatMult(gra,t(X)) #a%*%X
  ygam1=ygam2=ygam3=matrix(0,J,r)
  ygam2=t(grgamma[1,,])
  ygam3=t(grgamma[2,,])
  ygamat1=eigenMapMatMult(ygam1,t(X))
  ygamat2=eigenMapMatMult(ygam2,t(X))
  ygamat3=eigenMapMatMult(ygam3,t(X))
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
  aN=log(log(N))*log(2*J)
  BIC=-2*sum(lh)+l0norm*aN
  
  Bias=c(colSums(gra-Amat1)/10,colMeans(grd-Dmat1))
  RMSE=c(sqrt(colSums((gra-Amat1)^2)/10),sqrt(colMeans((grd-Dmat1)^2)))
  
  #output esimates and number of iterations
  return(list(est=cbind(gra,grd),Gamma=grgamma,mean1=Mu.gp1,mean2=Mu.gp2,mean3=Mu.gp3,Corr1=Sig.gp1,Corr2=Sig.gp2,Corr3=Sig.gp3,iter=iter,bic=BIC,bias=Bias,RMSE=RMSE))
}
#end of function

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

for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  r=2
  m=2
  lam=1
  eta.vec=seq(5,19,2)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    ptm <- proc.time()
    sim=ipest1(resp,m,r,eta,lam,eps =1e-3,max.tol=1e-7,NonUniform=T,gra00=gra00,grd00=grd00,grgamma00=grgamma00,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
    print(proc.time() - ptm)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Gammas[,,,k]=sim$Gamma
    biass[k,]=sim$bias
    RMSEs[k,]=sim$RMSE
  }
  
  kk=which.min(bics)
  eta.5[rep]=eta.vec[kk]
  Gammas.5[,,,rep]=Gammas[,,,kk]
  ADmat.5[,,rep]=ADmat[,,kk]
  biass.5[rep,]=biass[kk,]
  RMSEs.5[rep,]=RMSEs[kk,]
  print(ADmat.5[,,rep])
  print(eta.5[rep])
  print(Gammas.5[,,,rep])
  print(biass.5[rep,])
  print(RMSEs.5[rep,])
  write.csv(eta.5[rep],file = paste("eta6adapt_",rep))
  write.csv(ADmat.5[,,rep],file = paste("ADmat6adapt_",rep))
  write.csv(rbind(t(rbind(Gammas.5[c(1,2),1,3:11,rep])),t(rbind(Gammas.5[c(1,2),2,12:20,rep]))),file = paste("Gamma6adapt_",rep))
}


#########################################
#####                               #####
#####  Comparing with mirt LRT add  #####
#####                               #####
#########################################
N1=N2=N3=500 
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
Group01=c(rep('G1', N1), rep('G2', N2))
Group02=c(rep('G1', N1), rep('G3', N3))
N=N1+N2+N3
mirtCluster(4)

mirt.p.mat1=matrix(0,50,18) # the first 18 is the number of replications
bias.mirt1=matrix(0,50,3)
rmse.mirt1=matrix(0,50,3)
difrec.mirt.fn1=matrix(0,50,3) # DIF magnitude recovery with false negative included
mirt.p.mat2=matrix(0,50,18) # 18 is the number of replications
bias.mirt2=matrix(0,50,3)
rmse.mirt2=matrix(0,50,3)
difrec.mirt.fn2=matrix(0,50,3) # 
mirt.p.mat3=matrix(0,50,18) # 18 is the number of replications
bias.mirt3=matrix(0,50,3)
rmse.mirt3=matrix(0,50,3)
difrec.mirt.fn3=matrix(0,50,3) # 

for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #md.noncons0 <- multipleGroup(resp, cmodel, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  #omnibus
  # From Phil, when you test DIF on slope you should also test the intercept at the same time
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[1:r]))
  #md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:r]))
  mirt.p.mat1[(rep),1:9]=DIF(md.noncons0, which.par = c('a1'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat1[(rep),10:18]=DIF(md.noncons0, which.par = c('a2'), p.adjust = 'fdr',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-(which(mirt.p.mat1[rep,]<0.05)+2)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-c()]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[rep,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[c(4,5)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[c(12,13)]-0.5)))
  difrec.mirt.fn1[rep,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[c(4,5)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[c(12,13)]-1)))
  #difrec.mirt.fn1[rep,2]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G2$items[,1])[c(4,5,6,7,8,9)]-0.5),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G2$items[,2])[c(12,13,14,15,16,17)]-0.5)))
  #difrec.mirt.fn1[rep,3]= mean(c(abs((coef(md.refit0,simplify=T)$G1$items[,1]-coef(md.refit0,simplify=T)$G3$items[,1])[c(4,5,6,7,8,9)]-1),abs((coef(md.refit0,simplify=T)$G1$items[,2]-coef(md.refit0,simplify=T)$G3$items[,2])[c(12,13,14,15,16,17)]-1)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3]))
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[1:r]))
  mirt.p.mat2[(rep),1:9]=DIF(md.noncons01, which.par = c('a1'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat2[(rep),10:18]=DIF(md.noncons01, which.par = c('a2'), p.adjust = 'fdr',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  #md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  if ((sum(mirt.p.mat2[rep,]<0.05))==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp01)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(c(abs((coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[c(4,5)]-0.5),abs((coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[c(12,13)]-0.5)))
  #difrec.mirt.fn2[rep,2]= mean(c(abs((coef(md.refit01,simplify=T)$G1$items[,1]-coef(md.refit01,simplify=T)$G2$items[,1])[c(4,5,6,7,8,9)]-0.5),abs((coef(md.refit01,simplify=T)$G1$items[,2]-coef(md.refit01,simplify=T)$G2$items[,2])[c(12,13,14,15,16,17)]-0.5)))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp02)[1:r]))
  mirt.p.mat3[(rep),1:9]=DIF(md.noncons02, which.par = c('a1'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:11))[,"adj_pvals"]
  mirt.p.mat3[(rep),10:18]=DIF(md.noncons02, which.par = c('a2'), p.adjust = 'fdr',scheme = 'add',items2test=c(12:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','intercepts',colnames(resp02)[-(which(mirt.p.mat3[rep,]<0.05)+2)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]= mean(c(abs((coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[c(4,5)]-1),abs((coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[c(12,13)]-1)))
  #difrec.mirt.fn3[rep,3]= mean(c(abs((coef(md.refit02,simplify=T)$G1$items[,1]-coef(md.refit02,simplify=T)$G3$items[,1])[c(4,5,6,7,8,9)]-1),abs((coef(md.refit02,simplify=T)$G1$items[,2]-coef(md.refit02,simplify=T)$G3$items[,2])[c(12,13,14,15,16,17)]-1)))
  rep
}
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt1)
sqrt(colMeans(rmse.mirt1^2))
colMeans(difrec.mirt.fn1)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(18*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(18*14)
colMeans(bias.mirt2)
sqrt(colMeans(rmse.mirt2^2))
colMeans(difrec.mirt.fn2)

sum(mirt.p.mat3[1:24,c(2,3,10,11)]<0.05)/(24*4)
sum(mirt.p.mat3[1:24,-c(2,3,10,11)]<0.05)/(24*14)
colMeans(bias.mirt3)
sqrt(colMeans(rmse.mirt3^2))
colMeans(difrec.mirt.fn3)



# 60% DIF
sum(mirt.p.mat1[,c(2,3,4,5,6,7,10,11,12,13,14,15)]<0.05)/(50*12)
sum(mirt.p.mat1[,-c(2,3,4,5,6,7,10,11,12,13,14,15)]<0.05)/(50*8)
colMeans(bias.mirt1)
sqrt(colMeans(rmse.mirt1^2))
colMeans(difrec.mirt.fn1)

sum(mirt.p.mat2[,c(2,3,4,5,6,7,10,11,12,13,14,15)]<0.05)/(50*12)
sum(mirt.p.mat2[,-c(2,3,4,5,6,7,10,11,12,13,14,15)]<0.05)/(50*8)
colMeans(bias.mirt2)
sqrt(colMeans(rmse.mirt2^2))
colMeans(difrec.mirt.fn2)

sum(mirt.p.mat3[,c(2,3,4,5,6,7,10,11,12,13,14,15)]<0.05)/(50*12)
sum(mirt.p.mat3[,-c(2,3,4,5,6,7,10,11,12,13,14,15)]<0.05)/(50*8)
colMeans(bias.mirt3)
sqrt(colMeans(rmse.mirt3^2))
colMeans(difrec.mirt.fn3)

write.csv(mirt.p.mat1,file = "Sim7_LRTpvs1.csv")
write.csv(mirt.p.mat2,file = "Sim7_LRTpvs2.csv")
write.csv(mirt.p.mat3,file = "Sim7_LRTpvs3.csv")
write.csv(bias.mirt1,file = "Sim7_LRTbias1.csv")
write.csv(bias.mirt2,file = "Sim7_LRTbias2.csv")
write.csv(bias.mirt3,file = "Sim7_LRTbias3.csv")
write.csv(rmse.mirt1,file = "Sim7_LRTrmse1.csv")
write.csv(rmse.mirt2,file = "Sim7_LRTrmse2.csv")
write.csv(rmse.mirt3,file = "Sim7_LRTrmse3.csv")
write.csv(difrec.mirt.fn1,file = "Sim7_LRTfn1.csv")
write.csv(difrec.mirt.fn2,file = "Sim7_LRTfn2.csv")
write.csv(difrec.mirt.fn3,file = "Sim7_LRTfn3.csv")



