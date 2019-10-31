# confirmatory (test MC)
# confirmatory (test quadrature)
set.seed(11)
#data generation
L= 30 
m=4
r=2
N=1000
set.seed(11)
gra1.t=c(runif(1,1.1,2.8),0,runif((L-r)/2,1.1,2.8),rep(0,(L-r)/2))
gra2.t=c(0,runif(1,1.1,2.8),rep(0,(L-r)/2),runif((L-r)/2,1.1,2.8))
gra.t=cbind(gra1.t,gra2.t)

# we estimate d=a*b rather than the difficulty parameter b. in mirt p=1/1+exp(-(a*th+d))
#the difficulty parameter d is the same dimension as latent ability theta, but d is one-dimensional
grd1.t=runif(30,0.67,2)
grd2.t=runif(30,-0.67,0.67)
grd3.t=runif(30,-2,-0.67)
grd.t=cbind(grd1.t,grd2.t,grd3.t)
grip=cbind(gra.t,grd.t)

Sig.t=matrix(0.4,r,r)
diag(Sig.t)=1
theta=rmvnorm(N,mu=rep(0,r),Sigma = Sig.t)
mean(theta)
var(theta)
resp = simdata(as.matrix(gra.t),as.matrix(grd.t),itemtype = "graded",Theta=theta)
eps =1e-3
max.tol=1e-7
eta=0
G=100 
  #make sure the responses are coded from 1 instead of 0
  if(min(resp)==0)
  {
    resp <- resp+1
  }
  #sample size and test length
  N <- nrow(resp)
  J <- ncol(resp)
  ##### starting values for the first eta
  if (eta==0){
    if (m==2){
      grd <- matrix(0,J,1)
    }else
      grd <- matrix(rep(seq(1,-1,length.out = m-1),each=J),J)
    gra.res <- matrix(0,nrow = r,ncol = r)
    diag(gra.res)=1.5
    gra.ran=matrix(1,J-r,r)
    gra=rbind(gra.res,gra.ran)
    gra[3:16,2]=0
    gra[17:30,1]=0
    Sig=matrix(0.1,r,r)
    diag(Sig)=1
  } else {
    ##### starting values for after the first eta
    if (m==2){
      grd=matrix(es0$est2[,(r+1):(m+r-1)],J,1)
    } else{
      grd=es0$est[,(r+1):(m+r-1)]
    }
    gra=es0$est2[,1:r]
    Sig=es0$Corr
  }
  
  Pstar <- Qstar <- matrix(0,N*G,m-1) #array(0,dim=c(G,m-1,N))
  P<- matrix(0,N*G,m) #array(0,dim=c(G,m,N))
  df.a <- df.d <- df.Sig <- 1
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps & max(df.d)>eps & max(df.Sig)>eps)
  {
    aold <- gra
    dold <- grd
    Sig.old <- Sig
    X=rmvnorm(n = G*N, rep(0,r), Sig)
    
    #calculation of n_g 
    axmat=gra%*%t(X) #a%*%X
    pstar=array(0,dim = c(J,(m-1),G*N))
    p=array(0,dim = c(J,m,G*N))
    for (g in 1:(G*N))
    {
      pstar[,,g] = 1/(1+exp(-(grd+axmat[,g])))
      p[,,g] = t(-diff(rbind(rep(1,J),t(pstar[,,g]),rep(0,J))))
    }
    pij=matrix(double(J*G),J,G)
    LiA=matrix(double(N*G),N,G)
    
    for (i in 1:N)
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G)
        {
          pij[j,g] <- p[j,respi[j],(G*(i-1)+g)]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)
    }
    Pi = apply(LiA,1,sum)
    ng.mat= LiA/Pi
    
    #update mu hat and Sigma hat
    ng.vec=numeric(G*N)
    for (i in 1:N){
      ng.vec[((i-1)*G+1):(i*G)]=ng.mat[i,]
    }
    Sigg=matrix(0,r,r)
    for (g in 1:(G*N)){
      Sigg=Sigg+X[g,]%*%t(X[g,])*ng.vec[g]
    }
    Sig.hat=Sigg/N
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat))
    Tau.mat=matrix(rep(Tau,G*N),G*N,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    
    Sigg2=matrix(0,r,r)
    for (g in 1:(G*N)){
      Sigg2=Sigg2+Xstar[g,]%*%t(Xstar[g,])*ng.vec[g]
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
      rigk <- rLiA
      
      #M-step loop starts for item j
      miter <- 0
      add <- max.tol+1
      while(sum(abs(add))>max.tol)
      {
        miter <- miter+1
        for(g in 1:(N*G)){
          Pstar[g,] <- 1/(1+exp(-(d+rep(a*X[g,j]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        #calculating the score vector
        Dsco_g=numeric(m-1)
        for (i in 1:N){
          if (m==2){
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          } else {
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          }
        }
        Dsco=apply(Dsco_g,2,sum)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco_g=0
        for (i in 1:N){
          Asco_g = Asco_g+X[((i-1)*G+1):(i*G),j]*apply(rigk[i,,]*PQdif[((i-1)*G+1):(i*G),]/P[((i-1)*G+1):(i*G),],1,sum)
        }
        Asco=sum(Asco_g)
        minusgrad <- -c(Dsco,Asco)
        
        FI <- matrix(0,m,m)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng.vec*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng.vec*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          for (mm in 2:(m-1)){
            FI[mm,mm-1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
          }
          for (mm in 1:(m-2)){
            FI[mm,mm+1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
          }
        }
        
        FI[m,m] <- -sum(ng.vec*X[,j]^2*apply(PQdif^2/P,1,sum))
        
        
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
    ##item 3-16
    #calculation of r_jgk
    for (j in 3:16)
    {
      d <- grd[j,]
      a <- gra[j,1]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
      }
      rigk <- rLiA
      
      #M-step loop starts for item j
      miter <- 0
      add <- max.tol+1
      while(sum(abs(add))>max.tol)
      {
        miter <- miter+1
        for(g in 1:(N*G)){
          Pstar[g,] <- 1/(1+exp(-(d+rep(a*X[g,1]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        #calculating the score vector
        Dsco_g=numeric(m-1)
        for (i in 1:N){
          if (m==2){
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          } else {
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          }
        }
        Dsco=apply(Dsco_g,2,sum)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco_g=0
        for (i in 1:N){
          Asco_g = Asco_g+X[((i-1)*G+1):(i*G),1]*apply(rigk[i,,]*PQdif[((i-1)*G+1):(i*G),]/P[((i-1)*G+1):(i*G),],1,sum)
        }
        Asco=sum(Asco_g)
        minusgrad <- -c(Dsco,Asco)
        
        FI <- matrix(0,m,m)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng.vec*X[,1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng.vec*X[,1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          for (mm in 2:(m-1)){
            FI[mm,mm-1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
          }
          for (mm in 1:(m-2)){
            FI[mm,mm+1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
          }
        }
        
        FI[m,m] <- -sum(ng.vec*X[,1]^2*apply(PQdif^2/P,1,sum))
        
        
        add <- qr.solve(FI,minusgrad)
        
        d <- d+add[1:(m-1)]
        a <- a+add[m]
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[1]
      
      gra[j,1] <- a
      grd[j,] <- d
    }   
    ##item 17-30
    #calculation of r_jgk
    for (j in 17:30)
    {
      d <- grd[j,]
      a <- gra[j,2]
      rLiA <- array(double(N*G*m),dim = c(N,G,m))
      for(i in 1:N){
        for (g in 1:G){
          rLiA[i,g,resp[i,j]] <- LiA[i,g]
        }
        rLiA[i,,] <- rLiA[i,,]/Pi[i]
      }
      rigk <- rLiA
      
      #M-step loop starts for item j
      miter <- 0
      add <- max.tol+1
      while(sum(abs(add))>max.tol)
      {
        miter <- miter+1
        for(g in 1:(N*G)){
          Pstar[g,] <- 1/(1+exp(-(d+rep(a*X[g,2]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        #calculating the score vector
        Dsco_g=numeric(m-1)
        for (i in 1:N){
          if (m==2){
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          } else {
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          }
        }
        Dsco=apply(Dsco_g,2,sum)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco_g=0
        for (i in 1:N){
          Asco_g = Asco_g+X[((i-1)*G+1):(i*G),2]*apply(rigk[i,,]*PQdif[((i-1)*G+1):(i*G),]/P[((i-1)*G+1):(i*G),],1,sum)
        }
        Asco=sum(Asco_g)
        minusgrad <- -c(Dsco,Asco)
        
        FI <- matrix(0,m,m)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng.vec*X[,2]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            FI[m,mm] <- FI[mm,m] <- sum(ng.vec*X[,2]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
          }
          for (mm in 2:(m-1)){
            FI[mm,mm-1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
          }
          for (mm in 1:(m-2)){
            FI[mm,mm+1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
          }
        }
        
        FI[m,m] <- -sum(ng.vec*X[,2]^2*apply(PQdif^2/P,1,sum))
        
        
        add <- qr.solve(FI,minusgrad)
        
        d <- d+add[1:(m-1)]
        a <- a+add[m]
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau[2]
      
      gra[j,2] <- a
      grd[j,] <- d
    }   
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.Sig = abs(Sig.old-Sig)
    iter <- iter+1
  }
  
  
    colSums(gra-gra.t)/15
  sqrt(colSums((gra.t-gra)^2)/15)
  colSums(grd-grd.t)/30
  sqrt(colSums((grd.t-grd)^2)/30)
  Sig
  
  
  
  
  
  
  
  
  
  
  
  
