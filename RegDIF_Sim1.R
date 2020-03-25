library(MASS)
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("RESP1.csv",row.names = 1)

soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }
# Dataset #4 (2 Non-uniform DIF items per scale)
J=20

N1=N2=N3=500 
Group=c(rep('G1', N1), rep('G2', N2), rep('G3', N3))
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

#########################
##     For \eta>0      ##
#########################

########################
####### r=2, m=2 #######
########################

#########################
##  start of function  ##
#########################

ipest1 <- function(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=gra00,grd00=grd00,grbeta00=grbeta00,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
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
  grbeta=grbeta00
  Sig.gp1=Sig100; Sig.gp2=Sig200; Sig.gp3=Sig300
  grgamma=array(0,dim=c(r,r,J))
  
  #starting values for the MC samples form the prior
  Mu.gp1=mu100; Mu.gp2=mu200; Mu.gp3=mu300
  #X=theta #for check later
  
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  
  iter <- 0
  
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    #gammaold=grgamma
    betaold=grbeta
    A1=dmvnorm(X,Mu.gp1,Sig.gp1)
    A2=dmvnorm(X,Mu.gp2,Sig.gp2)
    A3=dmvnorm(X,Mu.gp3,Sig.gp3)
    
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
          pij[j,g] <- p1[j,as.numeric(respi[j]),g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A1
    }
    
    for (i in (N1+1):(N1+N2))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p2[j,as.numeric(respi[j]),g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A2
    }
    
    for (i in (N1+N2+1):(N1+N2+N3))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p3[j,as.numeric(respi[j]),g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A3
    }
    
    Pi = apply(LiA,1,sum)
    ng = apply(LiA/Pi,2,sum)
    ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
    ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
    ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
    
    #update mu hat and Sigma hat
    Mu.gp2=colSums(X*ng2)/N2
    Mu.gp3=colSums(X*ng3)/N3
    
    #update mu hat and Sigma hat
    Sigg1=Sigg2=Sigg3=matrix(0,r,r)
    for (g in 1:G){
      Sigg1=Sigg1+X[g,]%*%t(X[g,])*ng1[g]
    }
    Sig.hat1=Sigg1/N1
    
    for (g in 1:G){
      Sigg2=Sigg2+(X[g,]-Mu.gp2)%*%t(X[g,]-Mu.gp2)*ng2[g]
    }
    Sig.hat2=Sigg2/N2
    
    for (g in 1:G){
      Sigg3=Sigg3+(X[g,]-Mu.gp3)%*%t(X[g,]-Mu.gp3)*ng3[g]
    }
    Sig.hat3=Sigg3/N3
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat1))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    
    Sigg1=Sigg2=Sigg3=matrix(0,r,r)
    for (g in 1:G){
      Sigg1=Sigg1+Xstar[g,]%*%t(Xstar[g,])*ng1[g]
    }
    Sig.gp1=Sigg1/N1
    
    for (g in 1:G){
      Sigg2=Sigg2+(Xstar[g,]-Mu.gp2)%*%t(Xstar[g,]-Mu.gp2)*ng2[g]
    }
    Sig.gp2=Sigg2/N2
    
    for (g in 1:G){
      Sigg3=Sigg3+(Xstar[g,]-Mu.gp3)%*%t(Xstar[g,]-Mu.gp3)*ng3[g]
    }
    Sig.gp3=Sigg3/N3
    
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
      #gam=grgamma[,,j]
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
          Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y1%*%bet)))#+rep((y1%*%gam)%*%X[g,])
          P1[g,] <- -diff(c(1,Pstar1[g,],0))
          Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y2%*%bet)))#+rep((y2%*%gam)%*%X[g,])
          P2[g,] <- -diff(c(1,Pstar2[g,],0))
          Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y3%*%bet)))#+rep((y3%*%gam)%*%X[g,])
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
          Gamsco = c(apply(X[,1]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
        }
        if (m==2){
          Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        }else
          Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
          grad <- c(Dsco,Asco,Gamsco, Betsco)
        } else {
          minusgrad <- -c(Dsco,Asco, Betsco)
          grad <- c(Dsco,Asco, Betsco)
        }
        
        if (NonUniform==T){
          FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
          if (m==2){
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
              for (nn in m:(m+r-1)){
                FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
              }
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
          for (mm in m:(m+r-1)){
            for (nn in m:(m+r-1)){
              FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
            }
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
        bet0=bet+add[(m+r-1):(m+r)]
        for (mm in (m+r-1):(m+r)){
          add[mm]=soft(bet0[mm-m],-eta/FI[mm,mm])-bet[mm-m]
        }
        for (mm in (m+r-1):(m+r)){
          bet[mm-m]=soft(bet0[mm-m],-eta/FI[mm,mm])
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[1]
      
      gra[j,1] <- a
      grd[j,] <- d
      grbeta[j,] <- bet
    }
    
    # dimension2
    for (j in (r+J/2):J)
    {
      d <- grd[j,] 
      a <- gra[j,2]
      #gam=grgamma[,,j]
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
          Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y1%*%bet)))#+rep((y1%*%gam)%*%X[g,])
          P1[g,] <- -diff(c(1,Pstar1[g,],0))
          Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y2%*%bet)))#+rep((y2%*%gam)%*%X[g,])
          P2[g,] <- -diff(c(1,Pstar2[g,],0))
          Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y3%*%bet)))#+rep((y3%*%gam)%*%X[g,])
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
          Gamsco = c(apply(X[,2]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
        }
        if (m==2){
          Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
        }else
          Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
        if (NonUniform==T){
          minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
          grad <- c(Dsco,Asco,Gamsco, Betsco)
        } else {
          minusgrad <- -c(Dsco,Asco, Betsco)
          grad <- c(Dsco,Asco, Betsco)
        }
        
        if (NonUniform==T){
          FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
          if (m==2){
            for (mm in 1:(m-1)){
              FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
              for (nn in m:(m+r-1)){
                FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
              }
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
          for (mm in m:(m+r-1)){
            for (nn in m:(m+r-1)){
              FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
            }
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
        bet0=bet+add[(m+r-1):(m+r)]
        for (mm in (m+r-1):(m+r)){
          add[mm]=soft(bet0[mm-m],-eta/FI[mm,mm])-bet[mm-m]
        }
        for (mm in (m+r-1):(m+r)){
          bet[mm-m]=soft(bet0[mm-m],-eta/FI[mm,mm])
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[2]
      
      gra[j,2] <- a
      grd[j,] <- d
      grbeta[j,] <- bet
    }
    
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    iter <- iter+1
  }
  
  ###############
  # Refit model #
  ###############
  
  #get the sparsity structure
  sparsity=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  }
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps)
  {
    aold <- gra
    dold <- grd
    #gammaold=grgamma
    betaold=grbeta
    A1=dmvnorm(X,Mu.gp1,Sig.gp1)
    A2=dmvnorm(X,Mu.gp2,Sig.gp2)
    A3=dmvnorm(X,Mu.gp3,Sig.gp3)
    
    #calculation of n_g 
    axmat=gra%*%t(X) #a%*%X
    #ygam1=ygam2=ygam3=matrix(0,20,2)
    #for (j in 1:20){
    #  ygam1[j,]=y1%*%grgamma[,,j]
    #}
    #for (j in 1:20){
    #  ygam2[j,]=y2%*%grgamma[,,j]
    #}
    #for (j in 1:20){
    #  ygam3[j,]=y3%*%grgamma[,,j]
    #}
    #ygamat1=ygam1%*%t(X)
    #ygamat2=ygam2%*%t(X)
    #ygamat3=ygam3%*%t(X)
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
          pij[j,g] <- p1[j,as.numeric(respi[j]),g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A1
    }
    
    for (i in (N1+1):(N1+N2))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p2[j,as.numeric(respi[j]),g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A2
    }
    
    for (i in (N1+N2+1):(N1+N2+N3))
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p3[j,as.numeric(respi[j]),g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A3
    }
    
    Pi = apply(LiA,1,sum)
    ng = apply(LiA/Pi,2,sum)
    ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
    ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
    ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
    
    #update mu hat and Sigma hat
    Mu.gp2=colSums(X*ng2)/N2
    Mu.gp3=colSums(X*ng3)/N3
    
    #update mu hat and Sigma hat
    Sigg1=Sigg2=Sigg3=matrix(0,r,r)
    for (g in 1:G){
      Sigg1=Sigg1+X[g,]%*%t(X[g,])*ng1[g]
    }
    Sig.hat1=Sigg1/N1
    
    for (g in 1:G){
      Sigg2=Sigg2+(X[g,]-Mu.gp2)%*%t(X[g,]-Mu.gp2)*ng2[g]
    }
    Sig.hat2=Sigg2/N2
    
    for (g in 1:G){
      Sigg3=Sigg3+(X[g,]-Mu.gp3)%*%t(X[g,]-Mu.gp3)*ng3[g]
    }
    Sig.hat3=Sigg3/N3
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat1))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    
    Sigg1=Sigg2=Sigg3=matrix(0,r,r)
    for (g in 1:G){
      Sigg1=Sigg1+Xstar[g,]%*%t(Xstar[g,])*ng1[g]
    }
    Sig.gp1=Sigg1/N1
    
    for (g in 1:G){
      Sigg2=Sigg2+(Xstar[g,]-Mu.gp2)%*%t(Xstar[g,]-Mu.gp2)*ng2[g]
    }
    Sig.gp2=Sigg2/N2
    
    for (g in 1:G){
      Sigg3=Sigg3+(Xstar[g,]-Mu.gp3)%*%t(Xstar[g,]-Mu.gp3)*ng3[g]
    }
    Sig.gp3=Sigg3/N3
    
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
      #gam=grgamma[,,j]
      bet=grbeta[j,]
      l0normbetaj=0
      for(l in 1:2) 
      {
        l0normbetaj=l0normbetaj+(sparsity[j,l]!=0)
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
      if (l0normbetaj==0){
        miter <- 0
        add <- max.tol+1
        while(sum(abs(add))>max.tol)
        {
          miter <- miter+1
          for(g in 1:G){
            Pstar[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1]))))#+rep((y1%*%gam)%*%X[g,])
            P[g,] <- -diff(c(1,Pstar[g,],0))
          }
          Qstar <- 1-Pstar
          
          #calculating the score vector
          if (m==2){
            Dsco <- sum(apply(rgk/P,1,diff)*Pstar*Qstar)
          }else
            Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
          PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
          Asco = sum(X[,1]*apply(rgk*PQdif/P,1,sum))
          #  if (NonUniform==T){
          #  Gamsco = c(apply(X[,1]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
          #}
          
          if (NonUniform==T){
            minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
            grad <- c(Dsco,Asco,Gamsco, Betsco)
          } else {
            minusgrad <- -c(Dsco,Asco)
            grad <- c(Dsco,Asco)
          }
          
          if (NonUniform==T){
            FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
                for (nn in m:(m+r-1)){
                  FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
                }
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
            for (mm in m:(m+r-1)){
              for (nn in m:(m+r-1)){
                FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
              }
            }
            
          } else {
            FI <- matrix(0,m,m)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
                FI[m,mm] <- FI[mm,m] <- sum(ng*X[,1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
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
            
            FI[m,m] <- -sum(ng*X[,1]^2*apply(PQdif^2/P,1,sum))
          }
          
          
          add <- qr.solve(FI,minusgrad)
          d=d+add[1:(m-1)]
          a=a+add[m]
        }
        #end of M step loop
      } else {
        #M-step loop starts for item j
        miter <- 0
        add <- max.tol+1
        while(sum(abs(add))>max.tol)
        {
          miter <- miter+1
          for(g in 1:G){
            Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y1%*%bet)))#+rep((y1%*%gam)%*%X[g,])
            P1[g,] <- -diff(c(1,Pstar1[g,],0))
            Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y2%*%bet)))#+rep((y2%*%gam)%*%X[g,])
            P2[g,] <- -diff(c(1,Pstar2[g,],0))
            Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y3%*%bet)))#+rep((y3%*%gam)%*%X[g,])
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
          #if (NonUniform==T){
          #  Gamsco = c(apply(X[,1]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
          #}
          
          Betsco=numeric(l0normbetaj)
          if (m==2){
            if (l0normbetaj==1){
              Betsco[1] <- sparsity[j,]%*%c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
            } else {
              Betsco[1:l0normbetaj] <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
            }
            
          }else
            Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
          if (NonUniform==T){
            minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
            grad <- c(Dsco,Asco,Gamsco, Betsco)
          } else {
            minusgrad <- -c(Dsco,Asco, Betsco)
            grad <- c(Dsco,Asco, Betsco)
          }
          
          if (NonUniform==T){
            FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
                for (nn in m:(m+r-1)){
                  FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
                }
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
            for (mm in m:(m+r-1)){
              for (nn in m:(m+r-1)){
                FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
              }
            }
            
          } else {
            FI <- matrix(0,m+l0normbetaj,m+l0normbetaj)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
                FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
              }
              if (l0normbetaj==1){
                for (mm in (m+1):(m+l0normbetaj)){
                  FI[mm,mm] <- (sparsity[j,]%*%c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))
                  FI[m,mm] <- FI[mm,m] <- (sparsity[j,]%*%c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2]))))
                  FI[mm,m-1] <- FI[m-1,mm] <- (sparsity[j,]%*%c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))
                }
              } else {
                for (mm in (m+1):(m+l0normbetaj)){
                  FI[mm,mm] <- (c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))[mm-2]
                  FI[m,mm] <- FI[mm,m] <- (c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2]))))[mm-2]
                  FI[mm,m-1] <- FI[m-1,mm] <- (c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))[mm-2]
                }
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
          if (l0normbetaj==2){
            bet=bet+add[3:4]
          } else {
            if (sparsity[j,1]==1){
              bet=bet+c(add[(m+1):(m+l0normbetaj)],0)
            } else {
              bet=bet+c(0,add[(m+1):(m+l0normbetaj)])
            }
          }
        }
        #end of M step loop
        
      }
      
      #rescale a and d
      a=a*Tau[1]
      
      gra[j,1] <- a
      grd[j,] <- d
      grbeta[j,] <- bet
    }
    
    # dimension2
    for (j in (r+J/2):J)
    {
      d <- grd[j,] 
      a <- gra[j,2]
      #gam=grgamma[,,j]
      bet=grbeta[j,]
      l0normbetaj=0
      for(l in 1:2) 
      {
        l0normbetaj=l0normbetaj+(sparsity[j,l]!=0)
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
      if (l0normbetaj==0){
        miter <- 0
        add <- max.tol+1
        while(sum(abs(add))>max.tol)
        {
          miter <- miter+1
          for(g in 1:G){
            Pstar[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2]))))#+rep((y1%*%gam)%*%X[g,])
            P[g,] <- -diff(c(1,Pstar[g,],0))
          }
          Qstar <- 1-Pstar
          
          #calculating the score vector
          if (m==2){
            Dsco <- sum(apply(rgk/P,1,diff)*Pstar*Qstar)
          }else
            Dsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
          PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
          Asco = sum(X[,2]*apply(rgk*PQdif/P,1,sum))
          #  if (NonUniform==T){
          #  Gamsco = c(apply(X[,1]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
          #}
          
          if (NonUniform==T){
            minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
            grad <- c(Dsco,Asco,Gamsco, Betsco)
          } else {
            minusgrad <- -c(Dsco,Asco)
            grad <- c(Dsco,Asco)
          }
          
          if (NonUniform==T){
            FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
                for (nn in m:(m+r-1)){
                  FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
                }
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
            for (mm in m:(m+r-1)){
              for (nn in m:(m+r-1)){
                FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
              }
            }
            
          } else {
            FI <- matrix(0,m,m)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
                FI[m,mm] <- FI[mm,m] <- sum(ng*X[,2]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
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
            
            FI[m,m] <- -sum(ng*X[,2]^2*apply(PQdif^2/P,1,sum))
          }
          
          
          add <- qr.solve(FI,minusgrad)
          d=d+add[1:(m-1)]
          a=a+add[m]
        }
        #end of M step loop
      } else {
        miter <- 0
        add <- max.tol+1
        while(sum(abs(add))>max.tol)
        {
          miter <- miter+1
          for(g in 1:G){
            Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y1%*%bet)))#+rep((y1%*%gam)%*%X[g,])
            P1[g,] <- -diff(c(1,Pstar1[g,],0))
            Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y2%*%bet)))#+rep((y2%*%gam)%*%X[g,])
            P2[g,] <- -diff(c(1,Pstar2[g,],0))
            Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y3%*%bet)))#+rep((y3%*%gam)%*%X[g,])
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
          #if (NonUniform==T){
          #  Gamsco = c(apply(X[,2]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
          #}
          
          Betsco=numeric(l0normbetaj)
          if (m==2){
            if (l0normbetaj==1){
              Betsco[1] <- sparsity[j,]%*%c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
            } else {
              Betsco[1:l0normbetaj] <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
            }
          }else
            Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
          if (NonUniform==T){
            minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
            grad <- c(Dsco,Asco,Gamsco, Betsco)
          } else {
            minusgrad <- -c(Dsco,Asco, Betsco)
            grad <- c(Dsco,Asco, Betsco)
          }
          
          if (NonUniform==T){
            FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
                for (nn in m:(m+r-1)){
                  FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
                }
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
            for (mm in m:(m+r-1)){
              for (nn in m:(m+r-1)){
                FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
              }
            }
            
          } else {
            FI <- matrix(0,m+l0normbetaj,m+l0normbetaj)
            if (m==2){
              for (mm in 1:(m-1)){
                FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
                FI[m,mm] <- FI[mm,m] <- sum(ng1*X[,2]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,2]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,2]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
              }
              if (l0normbetaj==1){
                for (mm in (m+1):(m+l0normbetaj)){
                  FI[mm,mm] <- (sparsity[j,]%*%c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))
                  FI[m,mm] <- FI[mm,m] <- (sparsity[j,]%*%c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2]))))
                  FI[mm,m-1] <- FI[m-1,mm] <- (sparsity[j,]%*%c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))
                }
              } else {
                for (mm in (m+1):(m+l0normbetaj)){
                  FI[mm,mm] <- (c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))[mm-2]
                  FI[m,mm] <- FI[mm,m] <- (c(sum(ng2*X[,1]*Pstar2[,1]*Qstar2[,1]*(PQdif2[,1]/P2[,1]-PQdif2[,2]/P2[,2])),sum(ng3*X[,1]*Pstar3[,1]*Qstar3[,1]*(PQdif3[,1]/P3[,1]-PQdif3[,2]/P3[,2]))))[mm-2]
                  FI[mm,m-1] <- FI[m-1,mm] <- (c(-sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2])),-sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))))[mm-2]
                }
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
          if (l0normbetaj==2){
            bet=bet+add[3:4]
          } else {
            if (sparsity[j,1]==1){
              bet=bet+c(add[(m+1):(m+l0normbetaj)],0)
            } else {
              bet=bet+c(0,add[(m+1):(m+l0normbetaj)])
            }
          }
        }
        #end of M step loop
        
      }
      
      
      #rescale a and d
      a=a*Tau[2]
      
      gra[j,2] <- a
      grd[j,] <- d
      grbeta[j,] <- bet
    }
    
    
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    #df.Sig = abs(Sigold-Sig)
    iter <- iter+1
  }
  #output esimates and number of iterations
  
  ######################
  ###       BIC      ###
  ######################
  est.d=grd
  est.a=gra
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
  axmat=est.a%*%t(X) #a%*%X
  #ygam1=ygam2=ygam3=matrix(0,20,2)
  #for (j in 1:20){
  #  ygam1[j,]=y1%*%grgamma[,,j]
  #}
  #for (j in 1:20){
  #  ygam2[j,]=y2%*%grgamma[,,j]
  #}
  #for (j in 1:20){
  #  ygam3[j,]=y3%*%grgamma[,,j]
  #}
  #ygamat1=ygam1%*%t(X)
  #ygamat2=ygam2%*%t(X)
  #ygamat3=ygam3%*%t(X)
  grbeta1=t(y1%*%t(est.beta))
  grbeta2=t(y2%*%t(est.beta))
  grbeta3=t(y3%*%t(est.beta))
  pstar1=pstar2=pstar3=array(double(J*(m-1)*G),dim = c(J,(m-1),G))
  p1=p2=p3=array(double(J*m*G),dim = c(J,m,G))
  for (g in 1:G)
  {
    pstar1[,,g] = 1/(1+exp(-(est.d+axmat[,g]+grbeta1)))#+ygamat1[,g]
    p1[,,g] = t(-diff(rbind(rep(1,J),t(pstar1[,,g]),rep(0,J))))
  }
  for (g in 1:G)
  {
    pstar2[,,g] = 1/(1+exp(-(grd+axmat[,g]+grbeta2)))#+ygamat2[,g]
    p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
  }
  for (g in 1:G)
  {
    pstar3[,,g] = 1/(1+exp(-(grd+axmat[,g]+grbeta3)))#+ygamat3[,g]
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
        pij[j,g] <- p1[j,as.numeric(respi[j]),g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A1
  }
  
  for (i in (N1+1):(N1+N2))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p2[j,as.numeric(respi[j]),g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A2
  }
  
  for (i in (N1+N2+1):(N1+N2+N3))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p3[j,as.numeric(respi[j]),g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A3
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
      sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,])+ beta2))),0)))
    }
    sumoverk3=numeric(G)
    for(g in 1:G){
      sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,])+ beta3))),0)))
    }
    temp=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm=0
  for(i in 1:J) 
  {
    for(j in 1:2)
    {
      l0norm=l0norm+(est.beta[i,j]!=0)
    }
  }
  
  BIC=-2*sum(lh)+l0norm*log(N)

  Bias=c(colSums(gra-Amat1)/10,colMeans(grd-Dmat1))
  RMSE=c(sqrt(colSums((gra-Amat1)^2)/10),sqrt(colMeans((grd-Dmat1)^2)))
  
  return(list(est=cbind(gra,grd),Beta=grbeta,iter=iter,bic=BIC,bias=Bias,RMSE=RMSE,mean1=Mu.gp1,mean2=Mu.gp2,mean3=Mu.gp3,Corr1=Sig.gp1,Corr2=Sig.gp2,Corr3=Sig.gp3))
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
eta.1=numeric(reps)
#Gammas.1=array(double(2*J*m*50),dim = c(2,2,J,50))
Betas.1=array(double(J*2*reps),dim = c(J,2,reps))
ADmat.1=array(double(J*3*reps),dim = c(J,3,reps)) #a has 2 columns, d has 1 column
biass.1=matrix(0,reps,3)
RMSEs.1=matrix(0,reps,3)


#biass.lrt.1=matrix(0,50,3)
#RMSEs.lrt.1=matrix(0,50,3)
#Power.lrt.1=numeric(50)
#TypeI.lrt.1=numeric(50)

#biass.wald.1=matrix(0,50,3)
#RMSEs.wald.1=matrix(0,50,3)
#Power.wald.1=numeric(50)
#TypeI.wald.1=numeric(50)

###################################################################
#  starting values at eta=0 for later eta's for all replications  #
###################################################################

resp=responses[1:1500,]
s <- 'F1 = 1,3-11
          F2 = 2,12-20
          COV = F1*F2'
cmodel <- mirt.model(s)
md0 <- mirt(resp, cmodel, itemtype = "2PL")
gra0=coef(md0,simplify=T)$items[,1:r]
grd0=matrix(coef(md0,simplify=T)$items[,3],J,1)
grbeta0=matrix(0,J,2)
grbeta0[,1]=matrix(c(0,0,rep(0.7,18)),J,1)
grbeta0[,2]=matrix(c(0,0,rep(1.4,18)),J,1)
Sig0=matrix(c(1,coef(md0,simplify=T)$cov[2,1],coef(md0,simplify=T)$cov[2,1],1),2,2)


eta=0
eps =1e-3
max.tol=1e-7
NonUniform=F
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
grbeta=grbeta0
Sig.gp1=Sig.gp2=Sig.gp3=Sig0
grgamma=array(0,dim=c(r,r,J))

#starting values for the MC samples form the prior
Mu.gp1=Mu.gp2=Mu.gp3=rep(0,r)
#X=theta #for check later

Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1

iter <- 0

while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
{
  aold <- gra
  dold <- grd
  #gammaold=grgamma
  betaold=grbeta
  A1=dmvnorm(X,Mu.gp1,Sig.gp1)
  A2=dmvnorm(X,Mu.gp2,Sig.gp2)
  A3=dmvnorm(X,Mu.gp3,Sig.gp3)
  
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
        pij[j,g] <- p1[j,as.numeric(respi[j]),g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A1
  }
  
  for (i in (N1+1):(N1+N2))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p2[j,as.numeric(respi[j]),g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A2
  }
  
  for (i in (N1+N2+1):(N1+N2+N3))
  {
    respi=resp[i,]
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g] <- p3[j,as.numeric(respi[j]),g]
        
      }
    }
    LiA[i,]=apply(pij, 2, prod)*A3
  }
  
  Pi = apply(LiA,1,sum)
  ng = apply(LiA/Pi,2,sum)
  ng1=apply(LiA[1:N1,]/Pi[1:N1],2,sum)
  ng2=apply(LiA[(N1+1):(N1+N2),]/Pi[(N1+1):(N1+N2)],2,sum)
  ng3=apply(LiA[(N1+N2+1):(N1+N2+N3),]/Pi[(N1+N2+1):(N1+N2+N3)],2,sum)
  
  #update mu hat and Sigma hat
  Mu.gp2=colSums(X*ng2)/N2
  Mu.gp3=colSums(X*ng3)/N3
  
  #update mu hat and Sigma hat
  Sigg1=Sigg2=Sigg3=matrix(0,r,r)
  for (g in 1:G){
    Sigg1=Sigg1+X[g,]%*%t(X[g,])*ng1[g]
  }
  Sig.hat1=Sigg1/N1
  
  for (g in 1:G){
    Sigg2=Sigg2+(X[g,]-Mu.gp2)%*%t(X[g,]-Mu.gp2)*ng2[g]
  }
  Sig.hat2=Sigg2/N2
  
  for (g in 1:G){
    Sigg3=Sigg3+(X[g,]-Mu.gp3)%*%t(X[g,]-Mu.gp3)*ng3[g]
  }
  Sig.hat3=Sigg3/N3
  
  #scale 
  #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
  Tau=sqrt(diag(Sig.hat1))
  Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
  #q_g_star
  #X=(X-mu.hat.mat)/Tau.mat
  Xstar=X/Tau.mat
  
  Sigg1=Sigg2=Sigg3=matrix(0,r,r)
  for (g in 1:G){
    Sigg1=Sigg1+Xstar[g,]%*%t(Xstar[g,])*ng1[g]
  }
  Sig.gp1=Sigg1/N1
  
  for (g in 1:G){
    Sigg2=Sigg2+(Xstar[g,]-Mu.gp2)%*%t(Xstar[g,]-Mu.gp2)*ng2[g]
  }
  Sig.gp2=Sigg2/N2
  
  for (g in 1:G){
    Sigg3=Sigg3+(Xstar[g,]-Mu.gp3)%*%t(Xstar[g,]-Mu.gp3)*ng3[g]
  }
  Sig.gp3=Sigg3/N3
  
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
    #gam=grgamma[,,j]
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
        Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y1%*%bet)))#+rep((y1%*%gam)%*%X[g,])
        P1[g,] <- -diff(c(1,Pstar1[g,],0))
        Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y2%*%bet)))#+rep((y2%*%gam)%*%X[g,])
        P2[g,] <- -diff(c(1,Pstar2[g,],0))
        Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,1])+y3%*%bet)))#+rep((y3%*%gam)%*%X[g,])
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
        Gamsco = c(apply(X[,1]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
      }
      if (m==2){
        Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      }else
        Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      if (NonUniform==T){
        minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
        grad <- c(Dsco,Asco,Gamsco, Betsco)
      } else {
        minusgrad <- -c(Dsco,Asco, Betsco)
        grad <- c(Dsco,Asco, Betsco)
      }
      
      if (NonUniform==T){
        FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
            }
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
        for (mm in m:(m+r-1)){
          for (nn in m:(m+r-1)){
            FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
          }
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
      bet0=bet+add[(m+r-1):(m+r)]
      for (mm in (m+r-1):(m+r)){
        add[mm]=soft(bet0[mm-m],-eta/FI[mm,mm])-bet[mm-m]
      }
      for (mm in (m+r-1):(m+r)){
        bet[mm-m]=soft(bet0[mm-m],-eta/FI[mm,mm])
      }
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[1]
    
    gra[j,1] <- a
    grd[j,] <- d
    grbeta[j,] <- bet
  }
  
  # dimension2
  for (j in (r+J/2):J)
  {
    d <- grd[j,] 
    a <- gra[j,2]
    #gam=grgamma[,,j]
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
        Pstar1[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y1%*%bet)))#+rep((y1%*%gam)%*%X[g,])
        P1[g,] <- -diff(c(1,Pstar1[g,],0))
        Pstar2[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y2%*%bet)))#+rep((y2%*%gam)%*%X[g,])
        P2[g,] <- -diff(c(1,Pstar2[g,],0))
        Pstar3[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,2])+y3%*%bet)))#+rep((y3%*%gam)%*%X[g,])
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
        Gamsco = c(apply(X[,2]*apply(rgk2/P2*PQdif2,1,sum),2,sum),apply(X*apply(rgk3/P3*PQdif3,1,sum),2,sum))
      }
      if (m==2){
        Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      }else
        Betsco <- apply(t(apply(rgk/P,1,diff))*Pstar*Qstar,2,sum)
      if (NonUniform==T){
        minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
        grad <- c(Dsco,Asco,Gamsco, Betsco)
      } else {
        minusgrad <- -c(Dsco,Asco, Betsco)
        grad <- c(Dsco,Asco, Betsco)
      }
      
      if (NonUniform==T){
        FI <- matrix(0,m+r-1+r*2+(m-1)*2,m+r-1+r*2+(m-1)*2)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng1*(Pstar1[,mm]*Qstar1[,mm])^2*(1/P1[,mm]+1/P1[,mm+1])+ng2*(Pstar2[,mm]*Qstar2[,mm])^2*(1/P2[,mm]+1/P2[,mm+1])+ng3*(Pstar3[,mm]*Qstar3[,mm])^2*(1/P3[,mm]+1/P3[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng1*X[,nn-m+1]*Pstar1[,mm]*Qstar1[,mm]*(PQdif1[,mm]/P1[,mm]-PQdif1[,mm+1]/P1[,mm+1])+ng2*X[,nn-m+1]*Pstar2[,mm]*Qstar2[,mm]*(PQdif2[,mm]/P2[,mm]-PQdif2[,mm+1]/P2[,mm+1])+ng3*X[,nn-m+1]*Pstar3[,mm]*Qstar3[,mm]*(PQdif3[,mm]/P3[,mm]-PQdif3[,mm+1]/P3[,mm+1]))
            }
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
        for (mm in m:(m+r-1)){
          for (nn in m:(m+r-1)){
            FI[mm,nn] <- -sum(ng1*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif1^2/P1,1,sum)+ng2*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif2^2/P2,1,sum)+ng3*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif3^2/P3,1,sum))
          }
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
      bet0=bet+add[(m+r-1):(m+r)]
      for (mm in (m+r-1):(m+r)){
        add[mm]=soft(bet0[mm-m],-eta/FI[mm,mm])-bet[mm-m]
      }
      for (mm in (m+r-1):(m+r)){
        bet[mm-m]=soft(bet0[mm-m],-eta/FI[mm,mm])
      }
    }
    #end of M step loop
    
    #rescale a and d
    a=a*Tau[2]
    
    gra[j,2] <- a
    grd[j,] <- d
    grbeta[j,] <- bet
  }
  
  
  df.d <- abs(dold-grd)
  df.a <- abs(aold-gra)
  df.beta <- abs(betaold-grbeta)
  iter <- iter+1
}


gra00=gra
grd00=grd
grbeta00=grbeta
mu100= Mu.gp1
mu200= Mu.gp2
mu300= Mu.gp3
Sig100=Sig.gp1
Sig200=Sig.gp2
Sig300=Sig.gp3



for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  r=2
  m=2
  eta.vec=seq(21,48,3)
  bics=rep(0,length(eta.vec))
  ADmat=array(double(J*3*length(eta.vec)),dim = c(J,3,length(eta.vec)))
  #Gammas=array(double(2*J*m*length(eta.vec)),dim = c(2,2,J,length(eta.vec)))
  Betas=array(double(J*2*length(eta.vec)),dim = c(J,2,length(eta.vec)))
  biass=matrix(0,length(eta.vec),3)
  RMSEs=matrix(0,length(eta.vec),3)
  for (k in 1:length(eta.vec))
  {
    eta=eta.vec[k]
    sim=ipest1(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=gra00,grd00=grd00,grbeta00=grbeta00,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat=sim$est
    Betas[,,k]=sim$Beta
    biass[k,]=sim$bias
    RMSEs[k,]=sim$RMSE
  }
  
  kk=which.min(bics)
  
  eta.1[rep]=eta.vec[kk]
  #Gammas.13[,,,i]=Gammas[,,,kk]
  ADmat.1[,,rep]=ADmat[,,kk]
  Betas.1[,,rep]=Betas[,,kk]
  biass.1[rep,]=biass[kk,]
  RMSEs.1[rep,]=RMSEs[kk,]
  print(ADmat.1[,,rep])
  print(eta.1[rep])
  print(Betas.1[,,rep])
  print(biass.1[rep,])
  print(RMSEs.1[rep,])
}

sparsity.t=array(double(J*2*reps),dim = c(J,2,reps))
for (aa in 1:J){
  for (rr in 1:2){
    for (cc in 1:reps)
    sparsity.t[aa,rr,cc]=ifelse(Betas.1[aa,rr,cc]==0,0,1)
  }
}

power = c(sum(sparsity.t[c(4,5,12,13),1,])/(4*reps),sum(sparsity.t[c(4,5,12,13),2,])/(4*reps))
typeI = c(sum(sparsity.t[-c(4,5,12,13),1,])/(14*reps),sum(sparsity.t[-c(4,5,12,13),2,])/(14*reps))




