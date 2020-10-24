##################################################
##### Function: Reg_DIF                      #####
##################################################
##### Inputs
##### resp: response vector of length L, L > 1 
##### m: number of response categories
##### r: number of trait dimension 
##### y: number of examinee groups
##### N.vec: number of examinees in each group, 
##### vector of y, e.g. c(500,500,500)
##### mu.list: prior mean vector for each group, 
##### vector of r*y, e.g., the mean vector for 
##### three groups are (0,0), (0.1,0.02), (-0.05,0.03), 
##### then mu.list=c(0,0,0.1,0.02,-0.05,0.03)
##### Sig.list: prior covariance matrix for each group, 
##### matrix of r*y by r, each r by r matrix is the 
##### covariance matrix of each group
##### gra00: starting values of discrimination matrix, should include the information of known loading structure
##################################################                          
##### Outputs:
##### theta_hat: estimate of theta, vector of length p
##### std_error: standard error of the estimator, vector of length p 
##### detfi: determinant of observed Fisher Information, scalar                        
##################################################

#Reg_DIF <- function(resp,m,r,y,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=NULL,grd00=NULL,grbeta00=NULL,mu100=NULL,mu200=NULL,mu300=NULL,Sig100= NULL,Sig200= NULL,Sig300= NULL)
Reg_DIF <- function(resp,m,r,y,N.vec,loading.struc,eta,eps =1e-3,max.tol=1e-7,NonUniform=F,gra00=NULL,grd00=NULL,grbeta00=NULL,mu.list=NULL,Sig.list= NULL)
{
  #make sure the responses are coded from 1 instead of 0
  if(min(resp)==0)
  {
    resp <- resp+1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  # Gauss-Hermite quadrature nodes
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
  y.allgroup=rbind(rep(0,y-1),diag(y-1)) #y1, y2, y3
  # starting values
  gra=gra00
  grd=grd00
  grbeta=grbeta00 #matrix(0,J,2)
  grgamma=grgamma00 #array(0,dim=c((y-1),r,J))
  Sig.list=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.list=Mu.list #c(mu100,mu200,mu300)

  df.a <- df.d  <- df.gamma <- df.beta <- 1
  iter <- 0
  
  # regularied EM 
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps | max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta
    
    # E STEP
    LiA=E_step1(N.vec=c(500,500,500),X=X,y=y,y.allgroup=y.allgroup,Mu.list=Mu.list,Sig.list=Sig.list,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma)
    ng=ng.est(LiA=LiA,y=y,N.vec=c(500,500,500))
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.list
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-rep(Mu.est[((yy-1)*r+1):((yy-1)*r+r)])),((X-rep(Mu.est[((yy-1)*r+1):((yy-1)*r+r)]))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.list[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-rep(Mu.est[((yy-1)*r+1):((yy-1)*r+r)])),((Xstar-rep(Mu.est[((yy-1)*r+1):((yy-1)*r+r)]))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
   
    for (j in 1:J){
      rgk=rgk.est(j,Xijk,LiA,y,N.vec)
      Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
      P<- array(double(G*m*(y)),dim=c(G,m,y))
      M_step()
      #rescale a and d
      a=a*Tau[1]
      gra[j,1] <- a
      grd[j,] <- d
      grbeta[j,] <- bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    iter <- iter+1
  }
  # Re-estimation
  sparsity=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  }
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    #gammaold=grgamma
    betaold=grbeta
    
    # E STEP
    E_step1()
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
    
    for (j in 1:J){
      E_step2()
      M_step()
      #rescale a and d
      a=a*Tau[1]
      
      gra[j,1] <- a
      grd[j,] <- d
      grbeta[j,] <- bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    iter <- iter+1
  }
  # GIC
  E_step1()
  for (j in 1:J){
    E_step2()
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
  return(list(est=cbind(gra,grd),Beta=grbeta,iter=iter,bic=BIC,bias=Bias,RMSE=RMSE,mean1=Mu.gp1,mean2=Mu.gp2,mean3=Mu.gp3,Corr1=Sig.gp1,Corr2=Sig.gp2,Corr3=Sig.gp3))
}