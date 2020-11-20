##################################################
##### Function: Reg_DIF                      #####
##################################################
##### Inputs
##### resp: response vector of length L, L > 1 
##### m: number of response categories
##### r: number of trait dimension 
##### y: number of examinee groups y=3
##### N.vec: number of examinees in each group, 
##### vector of y, e.g. c(500,500,500)
##### mu.list: prior mean vector for each group, 
##### vector of r*y, e.g., the mean vector for 
##### three groups are (0,0), (0.1,0.02), (-0.05,0.03), 
##### then mu.list=c(0,0,0.1,0.02,-0.05,0.03)
##### Sig.list: prior covariance matrix for each group, 
##### matrix of r*y by r, each r by r matrix is the 
##### covariance matrix of each group
##### gra00: starting values of a, should include the information of known loading structure
##### grd00: starting values of d
##### grbeta00: starting values of beta
##### grgamma00: starting values of gamma
##################################################                          
##### Outputs:
##### theta_hat: estimate of theta, vector of length p
##### std_error: standard error of the estimator, vector of length p 
##### detfi: determinant of observed Fisher Information, scalar                        
##################################################

#Reg_DIF <- function(resp,m,r,y,N.vec=c(500,500,500),loading.struc,eta,eps =1e-3,max.tol=1e-7,NonUniform=T,gra00=NULL,grd00=NULL,grbeta00=NULL,mu.list=NULL,Sig.list= NULL)

NonUnif_Reg_DIF <- function(resp,m,r,y,N.vec,loading.struc,eta,eps =1e-3,max.tol=1e-7,gra00=NULL,grd00=NULL,grbeta00=NULL,mu.list=NULL,Sig.list= NULL)
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
    LiA=E_step1(N.vec=N.vec,X=X,y=y,y.allgroup=y.allgroup,Mu.list=Mu.list,Sig.list=Sig.list,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma)
    ng=ng.est(LiA=LiA,y=y,N.vec=N.vec)
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
      estj=M_step(j,grd,gra,grgamma,grbeta,max.tol,X,y.allgroup,y,G,eta)
      
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
  }
  # Re-estimation
  sparsity=grgamma
  for (j in 1:J){
    for (rr in 1:2){
      for (nn in 1:2){
        sparsity[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  
  gra=gra00
  grd=grd00
  grbeta=grbeta00
  grgamma=grgamma00*sparsity #array(0,dim=c((y-1),r,J))
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta
    
    # E STEP
    LiA=E_step1(N.vec=N.vec,X=X,y=y,y.allgroup=y.allgroup,Mu.list=Mu.list,Sig.list=Sig.list,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma)
    ng=ng.est(LiA=LiA,y=y,N.vec=N.vec)
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
      estj=M_step(j,grd,gra,grgamma,grbeta,max.tol,X,y.allgroup,y,G,eta=0)
      
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
  }
  # AIC BIC
  LiA=E_step1(N.vec=N.vec,X=X,y=y,y.allgroup=y.allgroup,Mu.list=Mu.list,Sig.list=Sig.list,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma)
  ng=ng.est(LiA=LiA,y=y,N.vec=N.vec)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgk.est(j,Xijk,LiA,y,N.vec)
    sumoverk1=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],gamjy=c(0,0),X=X)#G,rgky,aj,dj,gamjy,X
    sumoverk2=sumoverk(G=G,rgky=rgk[(2*G+1):(3*G),],aj=gra[j,],dj=grd[j,],gamjy=grgamma[1,,j],X=X)
    sumoverk3=sumoverk(G=G,rgky=rgk[(3*G+1):(4*G),],aj=gra[j,],dj=grd[j,],gamjy=grgamma[2,,j],X=X)
    temp=sum(sumoverk1)+sum(sumoverk2)+sum(sumoverk3)#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  for (j in 1:J){
    rgk=rgk.est(j,Xijk,LiA,y,N.vec)
    rgk1=rgk[(G+1):(2*G),]
    rgk2=rgk[(2*G+1):(3*G),]
    rgk3=rgk[(3*G+1):(4*G),]
    sumoverk1=numeric(G)
    for(g in 1:G){
      sumoverk1[g]=rgk1[g,]%*%log(-diff(c(1,1/(1+exp(-(grd[j,]+rep(gra[j,]%*%X[g,])))),0)))
    }
    sumoverk2=numeric(G)
    for(g in 1:G){
      sumoverk2[g]=rgk2[g,]%*%log(-diff(c(1,1/(1+exp(-(grd[j,]+rep(gra[j,]%*%X[g,])+rep(grgamma[1,,j]%*%X[g,])))),0)))
    }
    sumoverk3=numeric(G)
    for(g in 1:G){
      sumoverk3[g]=rgk3[g,]%*%log(-diff(c(1,1/(1+exp(-(grd[j,]+rep(gra[j,]%*%X[g,])+rep(grgamma[2,,j]%*%X[g,])))),0)))
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











sparsity=grbeta
for (j in 1:J){
  for (rr in 1:2){
    sparsity[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
  }
}