# fix thetas distribution and item parameters ad
library(MASS)
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
setwd('/Users/ruoyizhu/Documents/GitHub/RegDIF_SimData')
params=read.csv("Para2.csv",row.names = 1)
responses=read.csv("RESP2.csv",row.names = 1)

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

##############################################
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
theta.mirt=fscores(md00, method='MAP')
head(theta.mirt)

grbeta=grbeta000
gra=as.matrix(gra000)
rownames(gra) <- c()
grd=grd000

theta.ng=theta.mirt
theta.ng[which(theta.mirt[,1]<X1[1]),1]=-3
theta.ng[which(theta.mirt[,2]<X1[1]),2]=-3
for (gg in 1:(length(X1)-1)){
  theta.ng[which(theta.mirt[,1]>X1[gg]&theta.mirt[,1]<X1[gg+1]),1]=X1[gg]
  theta.ng[which(theta.mirt[,2]>X1[gg]&theta.mirt[,2]<X1[gg+1]),2]=X1[gg]
}
theta.ng[which(theta.mirt[,1]>X1[31]),1]=3
theta.ng[which(theta.mirt[,2]>X1[31]),2]=3
for (g in 1:G){
  ng[g]=sum(theta.ng[which(theta.ng[,1]==X[g,1]),2]==X[g,2])
}

theta.ng1=theta.mirt[1:500,]
theta.ng1[which(theta.mirt[1:500,1]<X1[1]),1]=-3
theta.ng1[which(theta.mirt[1:500,2]<X1[1]),2]=-3
for (gg in 1:(length(X1)-1)){
  theta.ng1[which(theta.mirt[1:500,1]>X1[gg]&theta.mirt[1:500,1]<X1[gg+1]),1]=X1[gg]
  theta.ng1[which(theta.mirt[1:500,2]>X1[gg]&theta.mirt[1:500,2]<X1[gg+1]),2]=X1[gg]
}
theta.ng1[which(theta.mirt[1:500,1]>X1[31]),1]=3
theta.ng1[which(theta.mirt[1:500,2]>X1[31]),2]=3
for (g in 1:G){
  ng1[g]=sum(theta.ng1[which(theta.ng1[,1]==X[g,1]),2]==X[g,2])
}

theta.ng2=theta.mirt[501:1000,]
theta.ng2[which(theta.mirt[501:1000,1]<X1[1]),1]=-3
theta.ng2[which(theta.mirt[501:1000,2]<X1[1]),2]=-3
for (gg in 1:(length(X1)-1)){
  theta.ng2[which(theta.mirt[501:1000,1]>X1[gg]&theta.mirt[501:1000,1]<X1[gg+1]),1]=X1[gg]
  theta.ng2[which(theta.mirt[501:1000,2]>X1[gg]&theta.mirt[501:1000,2]<X1[gg+1]),2]=X1[gg]
}
theta.ng2[which(theta.mirt[501:1000,1]>X1[31]),1]=3
theta.ng2[which(theta.mirt[501:1000,2]>X1[31]),2]=3
for (g in 1:G){
  ng2[g]=sum(theta.ng2[which(theta.ng2[,1]==X[g,1]),2]==X[g,2])
}

theta.ng3=theta.mirt[1001:1500,]
theta.ng3[which(theta.mirt[1001:1500,1]<X1[1]),1]=-3
theta.ng3[which(theta.mirt[1001:1500,2]<X1[1]),2]=-3
for (gg in 1:(length(X1)-1)){
  theta.ng3[which(theta.mirt[1001:1500,1]>X1[gg]&theta.mirt[1001:1500,1]<X1[gg+1]),1]=X1[gg]
  theta.ng3[which(theta.mirt[1001:1500,2]>X1[gg]&theta.mirt[1001:1500,2]<X1[gg+1]),2]=X1[gg]
}
theta.ng3[which(theta.mirt[1001:1500,1]>X1[31]),1]=3
theta.ng3[which(theta.mirt[1001:1500,2]>X1[31]),2]=3
for (g in 1:G){
  ng3[g]=sum(theta.ng3[which(theta.ng3[,1]==X[g,1]),2]==X[g,2])
}

for (g in 1:G){
  rgk[g,1]=sum(((resp[which(theta.ng[,1]==X[g,1]),j])[which(theta.ng[which(theta.ng[,1]==X[g,1]),2]==X[g,2])])==1)
  rgk[g,2]=sum(((resp[which(theta.ng[,1]==X[g,1]),j])[which(theta.ng[which(theta.ng[,1]==X[g,1]),2]==X[g,2])])==2)
}
resp1=resp[1:500,]
resp2=resp[501:1000,]
resp3=resp[1001:1500,]
for (g in 1:G){
  rgk1[g,1]=sum(((resp1[which(theta.ng1[,1]==X[g,1]),j])[which(theta.ng1[which(theta.ng1[,1]==X[g,1]),2]==X[g,2])])==1)
  rgk1[g,2]=sum(((resp1[which(theta.ng1[,1]==X[g,1]),j])[which(theta.ng1[which(theta.ng1[,1]==X[g,1]),2]==X[g,2])])==2)
}
for (g in 1:G){
  rgk2[g,1]=sum(((resp2[which(theta.ng2[,1]==X[g,1]),j])[which(theta.ng2[which(theta.ng2[,1]==X[g,1]),2]==X[g,2])])==1)
  rgk2[g,2]=sum(((resp2[which(theta.ng2[,1]==X[g,1]),j])[which(theta.ng2[which(theta.ng2[,1]==X[g,1]),2]==X[g,2])])==2)
}
for (g in 1:G){
  rgk3[g,1]=sum(((resp3[which(theta.ng3[,1]==X[g,1]),j])[which(theta.ng3[which(theta.ng3[,1]==X[g,1]),2]==X[g,2])])==1)
  rgk3[g,2]=sum(((resp3[which(theta.ng3[,1]==X[g,1]),j])[which(theta.ng3[which(theta.ng3[,1]==X[g,1]),2]==X[g,2])])==2)
}

######################################################################3


StartVals=read.csv("StartingValues2.csv",row.names = 1)
grbeta00=as.matrix(StartVals[,4:5])
rownames(grbeta00) <- c()
colnames(grbeta00) <- c()

# Fix a and d

eps =1e-3
max.tol=1e-7
 m=2
 r=2
 resp=responses[1:N,]
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
  grbeta=grbeta00
  grbeta=matrix(0,J,2)
  grbeta[3:20,1]=0.7
  grbeta[3:20,2]=1.4
 #####################
  gra=as.matrix(Amat1)
  rownames(gra) <- c()
  grd=Dmat1
  Mu.gp=c(0,0)
  Sig.gp=matrix(c(1,0.85,0.85,1),2,2)
  ####################
  gra=gra00
  grd=grd00
  Mu.gp1=mu100; Mu.gp2=mu200; Mu.gp3=mu300
  Sig.gp1=Sig100; Sig.gp2=Sig200; Sig.gp3=Sig300
  ####################
  
  
  Pstar <- Qstar <-Pstar1 <- Qstar1 <-Pstar2 <- Qstar2 <-Pstar3 <- Qstar3 <- matrix(double(G*(m-1)),G,m-1)
  P<-P1<-P2<-P3<- matrix(double(G*m),G,m)
  df.a <- df.d  <- df.beta <- 1
  
  iter <- 0
  
  while(max(df.beta)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    #gammaold=grgamma
    betaold=grbeta
    #A1=dmvnorm(X,Mu.gp,Sig.gp)
    #A2=dmvnorm(X,Mu.gp,Sig.gp)
    #A3=dmvnorm(X,Mu.gp,Sig.gp)
    A1=dmvnorm(X,Mu.gp1,Sig.gp1)
    A2=dmvnorm(X,Mu.gp2,Sig.gp2)
    A3=dmvnorm(X,Mu.gp3,Sig.gp3)
    
    #calculation of n_g 
    axmat=gra%*%t(X) #a%*%X
    grbeta1=t(y1%*%t(grbeta))
    grbeta2=t(y2%*%t(grbeta))
    grbeta3=t(y3%*%t(grbeta))
    pstar1=pstar2=pstar3=array(double(J*(m-1)*G),dim = c(J,(m-1),G))
    p1=p2=p3=array(double(J*m*G),dim = c(J,m,G))
    for (g in 1:G)
    {
      pstar1[,,g] = 1/(1+exp(-(grd+axmat[,g]+grbeta1)))
      p1[,,g] = t(-diff(rbind(rep(1,J),t(pstar1[,,g]),rep(0,J))))
    }
    for (g in 1:G)
    {
      pstar2[,,g] = 1/(1+exp(-(grd+axmat[,g]+grbeta2)))
      p2[,,g] = t(-diff(rbind(rep(1,J),t(pstar2[,,g]),rep(0,J))))
    }
    for (g in 1:G)
    {
      pstar3[,,g] = 1/(1+exp(-(grd+axmat[,g]+grbeta3)))
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
        
        Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      
        minusgrad <- -c(Betsco)
        grad <- c(Betsco)

        FI <- matrix(0,2,2)
        FI[1,1] <- -sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2]))
        FI[2,2] <- -sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))
    
        add <- qr.solve(FI,minusgrad)
        minusgrad[1]/FI[1,1]
        bet0=bet+add
        for (mm in 1:2){
          add[mm]=soft(bet0[mm],-eta/FI[mm,mm])-bet[mm]
        }
        for (mm in 1:2){
          bet[mm]=soft(bet0[mm],-eta/FI[mm,mm])
        }
        print(bet)
      }
      #end of M step loop
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
        Betsco <- c(sum(apply(rgk2/P2,1,diff)*Pstar2*Qstar2),sum(apply(rgk3/P3,1,diff)*Pstar3*Qstar3))
      
        minusgrad <- -c(Betsco)
  
          FI <- matrix(0,2,2)
          FI[1,1] <- -sum(ng2*(Pstar2[,1]*Qstar2[,1])^2*(1/P2[,1]+1/P2[,2]))
          FI[2,2] <- -sum(ng3*(Pstar3[,1]*Qstar3[,1])^2*(1/P3[,1]+1/P3[,2]))
              
        add <- qr.solve(FI,minusgrad)
        bet0=bet+add
        for (mm in 1:2){
          add[mm]=soft(bet0[mm],-eta/FI[mm,mm])-bet[mm]
        }
        for (mm in 1:2){
          bet[mm]=soft(bet0[mm],-eta/FI[mm,mm])
        }
      }
      print(bet)
      #end of M step loop

      grbeta[j,] <- bet
    }
    df.beta <- abs(betaold-grbeta)
    iter <- iter+1
  }
  
