#optim() function is very slow

library(MASS)
library(mirt) 
mirtCluster(4)
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
setwd('/Users/zhux0445/Documents/GitHub/RegDIF_SimData')
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para1.csv",row.names = 1)
responses=read.csv("RESP1lowcor.csv",row.names = 1)

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
Amat1=as.matrix(Amat1)
Amat2=as.matrix(Amat2)
Amat3=as.matrix(Amat3)
rownames(Amat1) <- c()
rownames(Amat2) <- c()
rownames(Amat3) <- c()
#test cor=0.25
#rep=1
#set.seed(rep*100)
#Theta=mvrnorm(n=N,mu=c(0,0),Sigma = matrix(c(1,0.25,0.25,1),2,2))
#Theta1=Theta[1:N1,]
#Theta2=Theta[(N1+1):(N1+N2),]
#Theta3=Theta[(N1+N2+1):(N1+N2+N3),]
#datasetR=simdata(Amat1,Dmat1,itemtype = "dich",Theta=Theta1)
#datasetF1=simdata(Amat2,Dmat2,itemtype = "dich",Theta=Theta2)
#datasetF2=simdata(Amat3,Dmat3,itemtype = "dich",Theta=Theta3)
#resp=rbind(datasetR,datasetF1,datasetF2)

#write.csv(cbind(gra00,grd00,grbeta00),file = "StartingValues1LowCor.csv")



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
      for(k in 1:m){
        rLiA[,,k]=Xijk[,j,k]*LiA/Pi
      }
      rgk <- apply(rLiA,c(2,3),sum)
      
      #for(i in 1:N){
      #  for (g in 1:G){
      #    rLiA[i,g,resp[i,j]] <- LiA[i,g]
      #  }
      #  rLiA[i,,] <- rLiA[i,,]/Pi[i]
      #}
      #rgk <- apply(rLiA,c(2,3),sum)
      
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
    
    for (j in 1:J)
    {
      d <- grd[j,] 
      a <- gra[j,]
      #gam=grgamma[,,j]
      bet=grbeta[j,]
      sparsityj=ifelse(c(d,a,bet),1,0)
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
        #loglikelihood function (for each item)
        sumoverk=numeric(G)
        obj <- function(x){
          for(g in 1:G){
            x1=(sparsityj*x)[1:(m-1)]  #x1 is d
            x2=(sparsityj*x)[m:(m+r-1)]  #x2 is a
            x3=(sparsityj*x)[m+r]  #x3 is beta1
            x4=(sparsityj*x)[m+r+1]  #x4 is beta2
            sumoverk[g]=rgk1[g,]%*%log(-diff(c(1,1/(1+exp(-(x1+rep(x2%*%X[g,])))),0)))+rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(x1+rep(x2%*%X[g,])+x3))),0)))+rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(x1+rep(x2%*%X[g,])+x4))),0)))
          }
          temp=sum(sumoverk)#sum over g
          return(-temp)
        }
        
        
        #solve for the parameters
        newp<-optim(c(d,a,bet), obj,  method = "L-BFGS-B")
        add<-newp$par-c(d,a,bet)
        d <- newp$par[1:(m-1)]
        a <- newp$par[m:(m+r-1)]
        bet=newp$par[(m+r):(m+r+1)]
        
        #print(miter)
        #plot(Pstar)
        #plot(PQdif)
        #print(minusgrad)
        #print(FI)
        #print(add)
        #print(c(d,a,bet))
      }
      #end of M step loop
      
      
      
      
      #rescale a and d
      a=a*Tau
      
      gra[j,] <- a
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
eta.2=numeric(reps)
#Gammas.1=array(double(2*J*m*50),dim = c(2,2,J,50))
Betas.2=array(double(J*2*reps),dim = c(J,2,reps))
ADmat.2=array(double(J*3*reps),dim = c(J,3,reps)) #a has 2 columns, d has 1 column
biass.2=matrix(0,reps,3)
RMSEs.2=matrix(0,reps,3)

#write.csv(cbind(gra00,grd00,grbeta00),file = "StartingValues2.csv")
StartVals=read.csv("StartingValues1LowCor.csv",row.names = 1)
gra00=as.matrix(StartVals[,1:2])
rownames(gra00) <- c()
grd00=matrix(StartVals[,3],20,1)
grbeta00=as.matrix(StartVals[,4:5])
rownames(grbeta00) <- c()
colnames(grbeta00) <- c()

# 6 dif per dim
mu100=c(0,0)
mu200=c(0.1472920, -0.0157996)
mu300=c(0.1100240, -0.1232052)
Sig100=matrix(c(1,0.2753316,0.2753316,1),2,2)
Sig200=matrix(c(1.3259608,0.3145355,0.3145355,1.1363796),2,2)
Sig300=matrix(c(1.2270710,0.2503095,0.2503095,1.0718629),2,2)

for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  r=2
  m=2
  eta.vec=seq(18,48,3)
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
    sim=ipest1(resp,m,r,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=gra00,grd00=grd00,grbeta00=grbeta00,mu100=mu100,mu200=mu200,mu300=mu300,Sig100=Sig100,Sig200=Sig200,Sig300=Sig300)
    bics[k]=sim$bic
    #Gammas[,,,k]=sim$Gamma
    ADmat[,,k]=sim$est
    Betas[,,k]=sim$Beta
    biass[k,]=sim$bias
    RMSEs[k,]=sim$RMSE
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
  write.csv(eta.2[rep],file = paste("eta1LowCor_",rep))
  write.csv(ADmat.2[,,rep],file = paste("ADmat1LowCor_",rep))
  write.csv(Betas.2[,,rep],file = paste("Beta1LowCor_",rep))
  write.csv(theta.dist[,,kk],file = paste("theta1LowCor_",rep))
}


#########################################
#####                               #####
#####  Comparing with mirt LRT add  #####
#####                               #####
#########################################

mirt.p.mat1=matrix(0,50,18) # 16 is the number of replications
bias.mirt1=matrix(0,50,3)
rmse.mirt1=matrix(0,50,3)
difrec.mirt.fn1=matrix(0,50,3) # DIF magnitude recovery with false negative included (only need the first column)
difrec.mirt.fn1.r=matrix(0,50,3) # DIF magnitude recovery for regularization method omnibus (only need the first column)
mirt.p.mat2=matrix(0,50,18) # 16 is the number of replications
bias.mirt2=matrix(0,50,3)
rmse.mirt2=matrix(0,50,3)
difrec.mirt.fn2=matrix(0,50,3) #(only need the second column)
mirt.p.mat3=matrix(0,50,18) # 16 is the number of replications
bias.mirt3=matrix(0,50,3)
rmse.mirt3=matrix(0,50,3)
difrec.mirt.fn3=matrix(0,50,3) #(only need the third column)

for (rep in 1:50){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  resp01=resp[1:(N1+N2),]
  resp02=rbind(resp[1:N1,],resp[(N1+N2+1):(N1+N2+N3),])
  s <- 'D1 = 1,3-11
          D2 = 2,12-20
          COV = D1*D2'
  #omnibus
  md.noncons0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:r]))
  mirt.p.mat1[(rep),]=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit0 <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-(which(mirt.p.mat1[rep,]<0.05)+2)]))
  #md.refit.r <- multipleGroup(resp, s, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[-c(4,5,12,13)]))
  bias.mirt1[rep,1:2]=colSums(coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt1[rep,1:2]=sqrt(colSums((coef(md.refit0,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt1[rep,3]=colMeans(coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt1[rep,3]=sqrt(colMeans((coef(md.refit0,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn1[rep,2:3]= c(mean(abs((coef(md.refit0,simplify=T)$G2$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,6,7,8,9,12,13,14,15,16,17)]-0.5)),mean(abs((coef(md.refit0,simplify=T)$G3$items[,3]-coef(md.refit0,simplify=T)$G1$items[,3])[c(4,5,6,7,8,9,12,13,14,15,16,17)]-1)))
  difrec.mirt.fn1[rep,1]=mean(c(difrec.mirt.fn1[rep,2],difrec.mirt.fn1[rep,3])) #omnibus DIF recovery in report (only need this)
  #difrec.mirt.fn1.r[rep,2:3]= c(mean(abs((coef(md.refit.r,simplify=T)$G2$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-0.5)),mean(abs((coef(md.refit.r,simplify=T)$G3$items[,3]-coef(md.refit.r,simplify=T)$G1$items[,3])[c(4,5,12,13)]-1)))
  #difrec.mirt.fn1.r[rep,1]=mean(c(difrec.mirt.fn1.r[rep,2],difrec.mirt.fn1.r[rep,3])) #omnibus DIF recovery in report (only need this)
  #ref vs focal1
  md.noncons01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:r]))
  mirt.p.mat2[(rep),]=DIF(md.noncons01, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  if ((sum(mirt.p.mat2[rep,]<0.05))==0){
    md.refit01 <-multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[1:J]))
  } else {
    md.refit01 <- multipleGroup(resp01, s, group = Group01,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp01)[-(which(mirt.p.mat2[rep,]<0.05)+2)]))
  }
  bias.mirt2[rep,1:2]=colSums(coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt2[rep,1:2]=sqrt(colSums((coef(md.refit01,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt2[rep,3]=colMeans(coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt2[rep,3]=sqrt(colMeans((coef(md.refit01,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn2[rep,2]= mean(abs((coef(md.refit01,simplify=T)$G2$items[,3]-coef(md.refit01,simplify=T)$G1$items[,3])[c(4,5,6,7,8,9,12,13,14,15,16,17)]-0.5))
  #ref vs focal2
  md.noncons02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[1:r]))
  mirt.p.mat3[(rep),]=DIF(md.noncons02, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(3:20))[,"adj_pvals"]
  md.refit02 <- multipleGroup(resp02, s, group = Group02,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp02)[-(which(mirt.p.mat3[rep,]<0.05)+2)]))
  bias.mirt3[rep,1:2]=colSums(coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)/10
  rmse.mirt3[rep,1:2]=sqrt(colSums((coef(md.refit02,simplify=T)$G1$items[,1:r]-Amat1)^2)/10)
  bias.mirt3[rep,3]=colMeans(coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)
  rmse.mirt3[rep,3]=sqrt(colMeans((coef(md.refit02,simplify=T)$G1$items[,3]-Dmat1)^2))
  difrec.mirt.fn3[rep,3]=mean(abs((coef(md.refit02,simplify=T)$G3$items[,3]-coef(md.refit02,simplify=T)$G1$items[,3])[c(4,5,6,7,8,9,12,13,14,15,16,17)]-1))
}
sum(mirt.p.mat1[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat1[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt1)
sqrt(colMeans((rmse.mirt1)^2))
colMeans(difrec.mirt.fn1)
colMeans(difrec.mirt.fn1.r)

sum(mirt.p.mat2[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat2[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt2)
sqrt(colMeans((rmse.mirt2)^2))
colMeans(difrec.mirt.fn2)

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

write.csv(mirt.p.mat1,file = "Sim4_LRTpvs1.csv")
write.csv(mirt.p.mat2,file = "Sim4_LRTpvs2.csv")
write.csv(mirt.p.mat3,file = "Sim4_LRTpvs3.csv")
write.csv(bias.mirt1,file = "Sim4_LRTbias1.csv")
write.csv(bias.mirt2,file = "Sim4_LRTbias2.csv")
write.csv(bias.mirt3,file = "Sim4_LRTbias3.csv")
write.csv(rmse.mirt1,file = "Sim4_LRTrmse1.csv")
write.csv(rmse.mirt2,file = "Sim4_LRTrmse2.csv")
write.csv(rmse.mirt3,file = "Sim4_LRTrmse3.csv")
write.csv(difrec.mirt.fn1,file = "Sim4_LRTfn1.csv")
write.csv(difrec.mirt.fn2,file = "Sim4_LRTfn2.csv")
write.csv(difrec.mirt.fn3,file = "Sim4_LRTfn3.csv")

sum(mirt.p.mat3[,c(2,3,10,11)]<0.05)/(50*4)
sum(mirt.p.mat3[,-c(2,3,10,11)]<0.05)/(50*14)
colMeans(bias.mirt3)
colMeans(rmse.mirt3)
colMeans(difrec.mirt.fn3)

