#MGRM MonteCarlo
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(irtoys)

soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

#start of function 
ipest2 <- function(resp,m,r,G,eta,eps =1e-3,max.tol=1e-7)
{
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
  
  #X=theta #for check later
  
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
    
    
    ##random part
    #calculation of r_jgk
    for (j in (r+1):J)
    {
      d <- grd[j,]
      a <- gra[j,]
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
          Pstar[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        
        #calculating the score vector
        Dsco_g=matrix(0,G,m-1)
        for (i in 1:N){
          if (m==2){
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          } else {
            Dsco_g <- Dsco_g+t(apply(rigk[i,,]/P[((i-1)*G+1):(i*G),],1,diff))*Pstar[((i-1)*G+1):(i*G),]*Qstar[((i-1)*G+1):(i*G),]
          }
        }
        Dsco=apply(Dsco_g,2,sum)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco_g=matrix(0,G,r)
        for (i in 1:N){
          Asco_g = Asco_g+X[((i-1)*G+1):(i*G),]*apply(rigk[i,,]*PQdif[((i-1)*G+1):(i*G),]/P[((i-1)*G+1):(i*G),],1,sum)
        }
        Asco=apply(Asco_g,2,sum)
        minusgrad <- -c(Dsco,Asco)
        grad <- c(Dsco,Asco)
        
        FI <- matrix(0,m+r-1,m+r-1)
        if (m==2){
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng.vec*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
            }
          }
        } else {
          for (mm in 1:(m-1)){
            FI[mm,mm] <- -sum(ng.vec*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
            for (nn in m:(m+r-1)){
              FI[nn,mm] <- FI[mm,nn] <- sum(ng.vec*X[,nn-m+1]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
            }
          }
          if (m>2){
            for (mm in 2:(m-1)){
              FI[mm,mm-1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
            }
            for (mm in 1:(m-2)){
              FI[mm,mm+1] <- sum(ng.vec*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
            }
          }
        }
        for (mm in m:(m+r-1)){
          for (nn in m:(m+r-1)){
            FI[mm,nn] <- -sum(ng.vec*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif^2/P,1,sum))
          }
        }
        add <- qr.solve(FI,minusgrad)
        d=d+add[1:(m-1)]
        a0=a+add[m:(m+r-1)]
        for (mm in m:(m+r-1)){
          add[mm]=soft(a0[mm-m+1],-eta/FI[mm,mm])-a[mm-m+1]
        }
        for (mm in m:(m+r-1)){
          a[mm-m+1]=soft(a0[mm-m+1],-eta/FI[mm,mm])
        }
      }
      #end of M step loop
      
      #rescale a and d
      a=a*Tau
      
      gra[j,] <- a
      grd[j,] <- d
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.Sig = abs(Sig.old-Sig)
    iter <- iter+1
  }
  #get the sparsity structure
  sparsity=gra
  for (j in 1:J){
    for (rr in 1:r){
      sparsity[j,rr]=ifelse(gra[j,rr]==0,0,1)
    }
  }
  colnamevec=c("F1","F2","F3","F4","F5")
  colnames(sparsity)=colnamevec[1:r]
  #Re-estimate by mirt
  COV <- matrix(TRUE, r,r)
  diag(COV)=FALSE
  cmodel=mirt.model(sparsity,COV=COV)
  md <- mirt(resp,cmodel,itemtype = "graded")
  estmirt=matrix(0,j,m+r-1)
  for(i in 1:j){
    estmirt[i,]=coef(md)[[i]]
  }
  gramirt=estmirt[,1:r]
  grdmirt=estmirt[,(r+1):(r+m-1)]
  #estimated latent trait
  theta.est=fscores(md,method = "EAP")
  #output esimates and number of iterations
  return(list(est=cbind(gra,grd),est2=cbind(gramirt,grdmirt),Corr=Sig,Corr12=coef(md)$GroupPars[4],iter=iter,score=theta.est))
}
#end of function

#start of simulation (Jiang et al.)
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

gra1.t=c(runif(1,1.1,2.8),0,0,runif(9,1.1,2.8),rep(0,18))
gra2.t=c(0,runif(1,1.1,2.8),0,rep(0,9),runif(9,1.1,2.8),rep(0,9))
gra3.t=c(0,0,runif(1,1.1,2.8),rep(0,18),runif(9,1.1,2.8))
gra.t=cbind(gra1.t,gra2.t,gra3.t)

# we estimate d=a*b rather than the difficulty parameter b. in mirt p=1/1+exp(-(a*th+d))
#the difficulty parameter d is the same dimension as latent ability theta, but d is one-dimensional.
grd.t=runif(L,-0.67,0.67)

grd1.t=runif(30,0.67,2)
grd2.t=runif(30,-0.67,0.67)
grd3.t=runif(30,-2,-0.67)
grd.t=cbind(grd1.t,grd2.t,grd3.t)
grip=cbind(gra.t,grd.t)
#(only adjacent boundary parameters with distance of at least 0.5 apart were retained)
#seed(1)
grip[4,6]=-1.7159522
grip[24,6]=-1.8159522


gra.t=grip[,1:3]
grd.t=grip[,4:6]

Sig.t=matrix(0.4,r,r)
diag(Sig.t)=1
theta=rmvnorm(N,rep(0,r),Sig.t)
mean(theta)
var(theta)
resp = simdata(as.matrix(gra.t),as.matrix(grd.t),itemtype = "dich",Theta=theta)
resp = simdata(as.matrix(gra.t),as.matrix(grd.t),itemtype = "graded",Theta=theta)

#item parameter estimation using the ipest function
set.seed(1)
es0=ipest2(resp,m,r,G=100,eta=0,eps =1e-3,max.tol=1e-7)
set.seed(1)
es5=ipest1(resp,m,r,eta=5,eps =1e-3,max.tol=1e-7)
set.seed(1)
es8=ipest1(resp,m,r,eta=8,eps =1e-3,max.tol=1e-7)
set.seed(1)
es9=ipest1(resp,m,r,eta=9,eps =1e-3,max.tol=1e-7)
set.seed(1)
es10=ipest2(resp,m,r,G=100,eta=10,eps =1e-3,max.tol=1e-7)
set.seed(1)
es15=ipest1(resp,m,r,eta=15,eps =1e-3,max.tol=1e-7)
set.seed(1)
es20=ipest1(resp,m,r,eta=20,eps =1e-3,max.tol=1e-7)
set.seed(1)
es22=ipest1(resp,m,r,eta=22,eps =1e-3,max.tol=1e-7)
set.seed(1)
es24=ipest1(resp,m,r,eta=24,eps =1e-3,max.tol=1e-7)
set.seed(1)
es25=ipest1(resp,m,r,eta=25,eps =1e-3,max.tol=1e-7)
set.seed(1)
es26=ipest1(resp,m,r,eta=26,eps =1e-3,max.tol=1e-7)
set.seed(1)
es30=ipest1(resp,m,r,eta=30,eps =1e-3,max.tol=1e-7)
set.seed(1)
es50=ipest1(resp,m,r,eta=50,eps =1e-3,max.tol=1e-7)
set.seed(1)
es70=ipest1(resp,m,r,eta=70,eps =1e-3,max.tol=1e-7)
es50$est




colMeans(es0$est-grip) 
sqrt(colMeans((grip-es0$est)^2))
colMeans(es0$est2-grip) 
sqrt(colMeans((grip-es0$est2)^2))
es50$Corr
es0$Corr12

set.seed(1)
es=ipest1(resp,m,r,eta=5,eps =1e-3,max.tol=1e-7)
colMeans(es0$est-grip) 
sqrt(colMeans((grip-es0$est)^2))
colMeans(es0$est2-grip) 
sqrt(colMeans((grip-es0$est2)^2))


est3.30 <- ipest1(resp,m,r,G=500,eta=30,eps =1e-3,max.tol=1e-7)
est3.0
colMeans(est3.0$est-grip)
sqrt(colMeans((grip-est3.0$est)^2))
est3.7 <- ipest1(resp,m,r,G=500,eta=7,eps =1e-3,max.tol=1e-7)
est3.7  
colMeans(est3.7$est-grip)
sqrt(colMeans((grip-est3.7$est)^2))
sum(est3.11$est==0)


######### BIC
est.d=matrix(es30$est2[,(r+1):(m+r-1)],L,m-1)
est.a=es30$est2[,1:r]
t=es30$score

if(min(resp)==0)
{
  resp <- resp+1
}

loglikelihood=numeric(N)
for (i in 1:N){
  pstar=1/(1+exp(-(matrix(rep(est.a%*%t[i,],m-1),J,byrow = F)+est.d)))
  p=-diff(rbind(rep(1,J),t(pstar),rep(0,J)))
  respi=resp[i,]
  pij=numeric(J)
  for (j in 1:J){
    pij[j]=p[respi[j],j]
  }
  loglikelihood[i]=log(prod(pij))
}
l0norm=0
for(i in 1:J) 
{
  for(j in 1:r)
  {
    l0norm=l0norm+(est.a[i,j]!=0)
  }
}

-2*sum(loglikelihood)+l0norm*log(N)



