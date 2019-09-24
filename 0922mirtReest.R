# Sep22 modification dichotomous case 
# Re-estimate by mirt
# Gauss-Hermite quadrature
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(irtoys)

mgauss.hermite <- function(n, mu, sigma, prune=NULL) {
  if(!all(dim(sigma) == length(mu)))stop("mu and sigma have nonconformable dimensions")
  dm  <- length(mu)
  gh  <-t(rbind(normal.qu(n)$quad.points,normal.qu(n)$quad.weights))
  #idx grows exponentially in n and dm
  idx <- as.matrix(expand.grid(rep(list(1:n),dm)))
  pts <- matrix(gh[idx,1],nrow(idx),dm)
  wts <- apply(matrix(gh[idx,2],nrow(idx),dm), 1, prod)
  ## prune
  if(!is.null(prune)) {
    qwt <- quantile(wts, probs=prune)
    pts <- pts[wts > qwt,]
    wts <- wts[wts > qwt]
  }
  ## rotate, scale, translate points
  eig <- eigen(sigma) 
  rot <- eig$vectors %*% diag(sqrt(eig$values))
  pts <- t(rot %*% t(pts) + mu)
  return(list(points=pts, weights=wts))
}


soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

########################
####### r=3, m=2 #######
########################

ipest1 <- function(resp,m,r,G=9,eta,eps =1e-3,max.tol=1e-7)
{
  #make sure the responses are coded from 1 instead of 0
  if(min(resp)==0)
  {
    resp <- resp+1
  }
  #sample size and test length
  N <- nrow(resp)
  J <- ncol(resp)
  
  #normal quadrature nodes X and weights A
  GH <- mgauss.hermite(G,rep(0,r),diag(r), prune=NULL)
  X <- GH$points
  A <- GH$weights
  
  ng <-  numeric(G^r)
  ##### starting values for the first eta
  #grd <- rep(0,30)
  #gra.res <- matrix(rep(0,r*r),nrow = r,ncol = r)
  #for(j in 1:r)
  #{
  #  gra.res[j,j]=1.5
  #}
  #gra.ran = matrix(c(rep(1.5,9),rep(0.5,18),rep(0.5,9),rep(1.5,9),rep(0.5,9),rep(0.5,18),rep(1.5,9)),nrow = J-r,ncol = r)
  #gra=rbind(gra.res,gra.ran)
  
  ##### starting values for after the first eta
  grd=est3.27.new$est[,4]
  gra=est3.27.new$est[,1:3]
  
  #grd=grd.t#for check E step
  #gra=gra.t
  
  Pstar <- Qstar <- matrix(double(G^r*(m-1)),G^r,m-1)
  P<- matrix(double(G^r*m),G^r,m)
  df.a <- df.d <- 1
  iter <- 0
  
  #start of the main loop
  while(max(df.a)>eps & max(df.d)>eps)
  {
    aold <- gra
    dold <- grd
    
    #calculation of n_g 
    axmat=gra%*%t(X) #a%*%X
    pstar=array(double(J*(m-1)*G^r),dim = c(J,(m-1),G^r))
    p=array(double(J*m*G^r),dim = c(J,m,G^r))
    for (g in 1:(G^r))
    {
      pstar[,,g] = 1/(1+exp(-(grd+axmat[,g])))
      p[,,g] = t(-diff(rbind(rep(1,30),t(pstar[,,g]),rep(0,30))))
    }
    pij=matrix(double(J*G^r),J,G^r)
    LiA=matrix(double(N*G^r),N,G^r)
    
    for (i in 1:N)
    {
      respi=resp[i,]
      for (j in 1:J)
      {
        for (g in 1:G^r)
        {
          pij[j,g] <- p[j,respi[j],g]
          
        }
      }
      LiA[i,]=apply(pij, 2, prod)*A
    }
    Pi = apply(LiA,1,sum)
    ng= apply(LiA/Pi,2,sum)
    
    ##Constraint part
    #calculation of r_jgk
    for (j in 1:r)
    {
      d <- grd[j]
      a <- gra[j,j]
      rLiA <- array(double(N*G^r*m),dim = c(N,G^r,m))
      for(i in 1:N){
        for (g in 1:G^r){
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
        for(g in 1:G^r){
          Pstar[g,] <- 1/(1+exp(-(d+rep(a*X[g,j]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        #calculating the score vector
        Dsco <- sum(apply(rgk/P,1,diff)*Pstar*Qstar)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco <- sum(X[,j]*apply(rgk*PQdif/P,1,sum))
        minusgrad <- -c(Dsco,Asco)
        
        FI <- matrix(0,m,m)
        for (mm in 1:(m-1)){
          FI[mm,mm] <- -sum(ng*(Pstar[,mm]*Qstar[,mm])^2*(1/P[,mm]+1/P[,mm+1]))
          FI[m,mm] <- FI[mm,m] <- sum(ng*X[,j]*Pstar[,mm]*Qstar[,mm]*(PQdif[,mm]/P[,mm]-PQdif[,mm+1]/P[,mm+1]))
        }
        if (m>2){
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
      
      gra[j,j] <- a
      grd[j] <- d
    }   
    
    
    ##random part
    #calculation of r_jgk
    for (j in (r+1):J)
    {
      d <- grd[j]
      a <- gra[j,]
      rLiA <- array(double(N*G^r*m),dim = c(N,G^r,m))
      for(i in 1:N){
        for (g in 1:G^r){
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
        for(g in 1:G^r){
          Pstar[g,] <- 1/(1+exp(-(d+rep(a%*%X[g,]))))
          P[g,] <- -diff(c(1,Pstar[g,],0))
        }
        Qstar <- 1-Pstar
        
        
        #calculating the score vector
        Dsco <- sum(apply(rgk/P,1,diff)*Pstar*Qstar)
        PQdif <- -t(apply(cbind(0,Pstar*Qstar,0),1, diff))
        Asco <- apply(X*apply(rgk/P*PQdif,1,sum),2,sum)
        minusgrad <- -c(Dsco,Asco)
        grad <- c(Dsco,Asco)
        
        FI <- matrix(0,m+r-1,m+r-1)
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
        for (mm in m:(m+r-1)){
          for (nn in m:(m+r-1)){
            FI[mm,nn] <- -sum(ng*X[,mm-m+1]*X[,nn-m+1]*apply(PQdif^2/P,1,sum))
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
      
      gra[j,] <- a
      grd[j] <- d
    }   
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    iter <- iter+1
  }
  
  #get the sparsity structure
  sparsity=gra
  for (j in 1:J){
    for (rr in 1:r){
      sparsity[j,rr]=ifelse(gra[j,rr]==0,0,1)
    }
  }
  colnames(sparsity)=c("F1","F2","F3")
  #Re-estimate by mirt
  cmodel=mirt.model(sparsity)
  md <- mirt(resp,cmodel,itemtype = "graded")
  estmirt=matrix(0,j,m+r-1)
  for(i in 1:j){
    estmirt[i,]=coef(md)[[i]]
  }
  gramirt=estmirt[,1:3]
  grdmirt=estmirt[,4]
 
  #output esimates and number of iterations
  return(list(est=cbind(gra,grd),est2=cbind(gramirt,grdmirt),iter=iter))
}
#end of function

#start of simulation (Jiang et al.)
set.seed(1)
#data generation
L= 30 
m=2
r=3
N=1000
set.seed(1)
gra1.t=c(runif(1,1.1,2.8),0,0,runif(9,1.1,2.8),rep(0,18))
gra2.t=c(0,runif(1,1.1,2.8),0,rep(0,9),runif(9,1.1,2.8),rep(0,9))
gra3.t=c(0,0,runif(1,1.1,2.8),rep(0,18),runif(9,1.1,2.8))
gra.t=cbind(gra1.t,gra2.t,gra3.t)

# we estimate d=a*b rather than the difficulty parameter b. 
#the difficulty parameter d is the same dimension as latent ability theta, but d is one-dimensional.
grd1.t=rep(0,30)
grd.t=grd1.t
grip=cbind(gra.t,grd.t)

#Sig=matrix(c(1,0.2,0.2,0.2,1,0.2,0.2,0.2,1),3,3)
Sig=diag(r)
theta=rmvnorm(N,mu=rep(0,r),Sigma = Sig)
resp = simdata(as.matrix(gra.t),as.matrix(grd.t),itemtype = "dich",Theta=theta)

#item parameter estimation using the ipest function
set.seed(1)
est3.28.new <- ipest1(resp,m,r,G=9,eta=28,eps =1e-3,max.tol=1e-7)
est3.28.new
apply(est3.28.new$est2-grip,2,sum)[1:3]/10
sqrt(apply((grip-est3.28.new$est2)^2,2,sum)[1:3]/10)
set.seed(1)
est3.50 <- ipest1(resp,m,r,G=9,eta=50,eps =1e-3,max.tol=1e-7)
est3.10
colMeans(gra-grip[,1:3])
sqrt(colMeans((grip-est3.1$est)^2))
sqrt(colMeans((grip[,1:3]-gra)^2))
sum(est3.11$est==0)




######### BIC

est.d=est3.28.new$est2[,4]
est.a=est3.28.new$est2[,1:3]


if(min(resp)==0)
{
  resp <- resp+1
}
G=9
r=3
GH <- mgauss.hermite(G,rep(0,r),diag(r), prune=NULL)
X <- GH$points
A <- GH$weights
J=L

##compute ng
axmat=est.a%*%t(X) #a%*%X
pstar=array(double(J*(m-1)*G^r),dim = c(J,(m-1),G^r))
p=array(double(J*m*G^r),dim = c(J,m,G^r))
for (g in 1:G^r)
{
  pstar[,,g] = 1/(1+exp(-(est.d+axmat[,g])))
  p[,,g] = t(-diff(rbind(rep(1,30),t(pstar[,,g]),rep(0,30))))
}
pij=matrix(double(J*G^r),J,G^r)
LiA=matrix(double(N*G^r),N,G^r)


for (i in 1:N)
{
  respi=resp[i,]
  for (j in 1:J)
  {
    for (g in 1:G^r)
    {
      pij[j,g] <- p[j,respi[j],g]
      
    }
  }
  LiA[i,]=apply(pij, 2, prod)*A
}
Pi = apply(LiA,1,sum)
ng= apply(LiA/Pi,2,sum)


lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
for(j in 1:J)
{
  d <- est.d[j]
  a <- est.a[j,]
  rLiA <- array(double(N*G^r*m),dim = c(N,G^r,m))
  for(i in 1:N){
    for (g in 1:G^r){
      rLiA[i,g,resp[i,j]] <- LiA[i,g]
    }
    rLiA[i,,] <- rLiA[i,,]/Pi[i]
  }
  rgk <- apply(rLiA,c(2,3),sum)
  
  
  sumoverk=numeric(G^r)
  for(g in 1:G^r){
    sumoverk[g]=rgk[g,]%*%log(-diff(c(1,1/(1+exp(-(d+rep(a%*%X[g,])))),0)))
  }
  temp=sum(sumoverk)#-eta*norm(as.matrix(x2), type = "1")##sum over g
  lh[j]=temp
}
l0norm=0
for(i in 1:J) 
{
  for(j in 1:r)
  {
    l0norm=l0norm+(est.a[i,j]!=0)
  }
}

-2*sum(lh)+l0norm*log(N)

