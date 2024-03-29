##### Function: Reg_DIF                      #####
##################################################
##### Inputs
##### resp: response vector of length L, L > 1 
##### m: number of response categories
##### r: number of trait dimension 
##### y: number of examinee groups y=3
##### N.vec: number of examinees in each group, 
##### vector of y, e.g. N.vec=c(500,500,500)
##### Mu.list: prior mean vector for each group, 
##### vector of r*y, e.g., the mean vector for 
##### three groups are (0,0), (0.1,0.02), (-0.05,0.03), 
##### then Mu.list=c(0,0,0.1,0.02,-0.05,0.03);Mu.list=c(mu100,mu200,mu300)
##### Sig.list: prior covariance matrix for each group, 
##### matrix of r*y by r, each r by r matrix is the 
##### covariance matrix of each group e.g. Sig100=matrix(c(1,0.8452613,0.8452613,1),2,2)
#####                                      Sig200=matrix(c(1.179328,1.065364,1.065364,1.179328),2,2);Sig300=matrix(c(0.9202015,0.8908855,0.8908855,0.9202015),2,2)
#####                                      Sig.list=rbind(Sig100,Sig200,Sig300)
##### gra00: starting values of a, should include the information of known loading structure
##### grd00: starting values of d
##### grbeta00: starting values of beta; Items with starting value zero are used as anchor. grbeta00=matrix(0,J,2)
##### grgamma00: starting values of gamma. Items with starting value zero are used as anchor.
##################################################                          
##### Outputs:
##### est: estimated parameter a and d.
##### Gamma: estimated parameter gamma
##### Beta: estimated parameter beta
##### iter: number of EM cycles
##### bic: BIC
##### means: estimated mean vector for each group
##### Covs: estimated covariance matrix for each group
##################################################
# 4

Reg_DIF <- function(resp,m,r,y,N.vec,eta,eps =1e-3,max.tol=1e-7,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
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
  Sig.est=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.est=Mu.list #c(mu100,mu200,mu300)
  
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
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
     #LiA=E_step0(resp=resp,N.vec=N.vec,X=X,y=y,G=G,y.allgroup=y.allgroup,Mu.list=c(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6]),Sig.list=rbind(Sig.est[1:2,],Sig.est[3:4,],Sig.est[5:6,]),gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta)
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r,J=J)
      #gra[j,] <- estj$a*Tau  # re-scale a and gamma
      #grd[j,] <- estj$d
      #grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # Re-estimation
  sparsity1=grgamma
  for (j in 1:J){
    for (rr in 1:2){
      for (nn in 1:2){
        sparsity1[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  } 
  
  gra=gra00
  grd=grd00
  grgamma=grgamma00*sparsity1
  grbeta=grbeta00*sparsity2

  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta
    
    # E STEP
    Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
    Sig.est.slice=array(0,c(2,2,3))
    Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
    LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
      P<- array(double(G*m*(y)),dim=c(G,m,y))
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0)
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r,J=J)
      #gra[j,] <- estj$a*Tau  # re-scale a and gamma
      #grd[j,] <- estj$d
      #grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(2, iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
  Sig.est.slice=array(0,c(2,2,3))
  Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
  LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk1=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=c(0,0),X=X)#G,rgky,aj,dj,gamjy,X
    sumoverk2=sumoverk(G=G,rgky=rgk[(2*G+1):(3*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,1],gamjy=grgamma[1,,j],X=X)
    sumoverk3=sumoverk(G=G,rgky=rgk[(3*G+1):(4*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,2],gamjy=grgamma[2,,j],X=X)
    temp=sumoverk1+sumoverk2+sumoverk3#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:2){
      for(k in 1:2){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:2){
      l0norm2=l0norm2+(grbeta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2
  
  BIC=-2*sum(lh)+l0norm*log(N)
  #Mu.gp1=Mu.est[1:2];Mu.gp2=Mu.est[3:4];Mu.gp3=Mu.est[5:6]
  #Sig.gp1=Sig.est[1:2,];Sig.gp2=Sig.est[3:4,];Sig.gp3=Sig.est[5:6,]
  return(list(est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta,iter=iter,bic=BIC, means=Mu.est,Covs=Sig.est))
}

# 5

Reg_EMM_DIF <- function(resp,m,r,y,N.vec,eta,eps =1e-3,max.tol=1e-7,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
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
  Sig.est=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.est=Mu.list #c(mu100,mu200,mu300)
  
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
    Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
    Sig.est.slice=array(0,c(2,2,3))
    Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
    #LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta)
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r,J=J)
      #gra[j,] <- estj$a*Tau  # re-scale a and gamma
      #grd[j,] <- estj$d
      #grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj$bet
    }
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      #estj2=M_step2(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta)
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0)
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r,J=J)
      #gra[j,] <- estj$a*Tau  # re-scale a and gamma
      #grd[j,] <- estj$d
      #grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
  Sig.est.slice=array(0,c(2,2,3))
  Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
  LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk1=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=c(0,0),X=X)#G,rgky,aj,dj,gamjy,X
    sumoverk2=sumoverk(G=G,rgky=rgk[(2*G+1):(3*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,1],gamjy=grgamma[1,,j],X=X)
    sumoverk3=sumoverk(G=G,rgky=rgk[(3*G+1):(4*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,2],gamjy=grgamma[2,,j],X=X)
    temp=sumoverk1+sumoverk2+sumoverk3#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:2){
      for(k in 1:2){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:2){
      l0norm2=l0norm2+(grbeta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2
  
  BIC=-2*sum(lh)+l0norm*log(N)
  #Mu.gp1=Mu.est[1:2];Mu.gp2=Mu.est[3:4];Mu.gp3=Mu.est[5:6]
  #Sig.gp1=Sig.est[1:2,];Sig.gp2=Sig.est[3:4,];Sig.gp3=Sig.est[5:6,]
  return(list(est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta,iter=iter,bic=BIC, means=Mu.est,Covs=Sig.est))
}


# 6

Reg_Adaptive_DIF <- function(resp,m,r,y,N.vec,eta,lam,eps =1e-3,max.tol=1e-7,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
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
  Sig.est=Sig.list #rbind(Sig100,Sig200,Sig300)
  Mu.est=Mu.list #c(mu100,mu200,mu300)
  
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
    Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
    Sig.est.slice=array(0,c(2,2,3))
    Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
    #LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=lam)
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      #estj=M_step_Adapt(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=1,r=r,J=J)
      #gra[j,] <- estj$a*Tau  # re-scale a and gamma
      #grd[j,] <- estj$d
      #grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # Re-estimation
  sparsity1=grgamma
  for (j in 1:J){
    for (rr in 1:2){
      for (nn in 1:2){
        sparsity1[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  } 
  
  gra=gra00
  grd=grd00
  grgamma=grgamma00*sparsity1
  grbeta=grbeta00*sparsity2
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
  {
    aold <- gra
    dold <- grd
    gammaold=grgamma
    betaold=grbeta
    
    # E STEP
    Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
    Sig.est.slice=array(0,c(2,2,3))
    Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
    #LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
    ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
    #update mu hat and Sigma hat
    Mu.est=numeric(r*y)
    for (yy in 2:y){
      Mu.est[((yy-1)*r+1):((yy-1)*r+r)]=colSums(X*ng[(yy*G+1):(yy*G+G)])/N.vec[yy]
    }
    #update Sigma hat
    Sig.hat.allgrp=Sig.est
    for (yy in 1:y){
      Sig.hat.allgrp[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((X-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    #scale 
    #mu.hat.mat=matrix(rep(mu.hat,G),G,r,byrow = T)
    Tau=sqrt(diag(Sig.hat.allgrp[1:r,]))
    Tau.mat=matrix(rep(Tau,G),G,r,byrow = T)
    #q_g_star
    #X=(X-mu.hat.mat)/Tau.mat
    Xstar=X/Tau.mat
    for (yy in 1:y){
      Sig.est[((yy-1)*r+1):((yy-1)*r+r),]=eigenMapMatMult(t(Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G))),((Xstar-cbind(rep(Mu.est[((yy-1)*r+1)],G),rep(Mu.est[((yy-1)*r+r)],G)))*ng[(yy*G+1):(yy*G+G)]))/N.vec[yy]
    }
    
    
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
      P<- array(double(G*m*(y)),dim=c(G,m,y))
      estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=lam)
      gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      grd[j,] <- estj[1:(m-1)]
      grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      #estj=M_step_Adapt(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=1,r=r,J=J)
      #gra[j,] <- estj$a*Tau  # re-scale a and gamma
      #grd[j,] <- estj$d
      #grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(2,iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=rbind(Mu.est[1:2],Mu.est[3:4],Mu.est[5:6])
  Sig.est.slice=array(0,c(2,2,3))
  Sig.est.slice[,,1]=Sig.est[1:2,];Sig.est.slice[,,2]=Sig.est[3:4,];Sig.est.slice[,,3]=Sig.est[5:6,]
  #LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
  LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk1=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=c(0,0),X=X)#G,rgky,aj,dj,gamjy,X
    sumoverk2=sumoverk(G=G,rgky=rgk[(2*G+1):(3*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,1],gamjy=grgamma[1,,j],X=X)
    sumoverk3=sumoverk(G=G,rgky=rgk[(3*G+1):(4*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,2],gamjy=grgamma[2,,j],X=X)
    temp=sumoverk1+sumoverk2+sumoverk3#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:2){
      for(k in 1:2){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:2){
      l0norm2=l0norm2+(grbeta[i,j]!=0)
    }
  }
  l0norm=l0norm1+l0norm2
  
  BIC=-2*sum(lh)+l0norm*log(N)
  #Mu.gp1=Mu.est[1:2];Mu.gp2=Mu.est[3:4];Mu.gp3=Mu.est[5:6]
  #Sig.gp1=Sig.est[1:2,];Sig.gp2=Sig.est[3:4,];Sig.gp3=Sig.est[5:6,]
  return(list(est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta,iter=iter,bic=BIC, means=Mu.est,Covs=Sig.est))
}
