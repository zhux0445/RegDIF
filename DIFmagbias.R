library(GLDEX)
reps=50
Gamma.8=array(double(18*2*reps),dim = c(18,2,reps))
power=matrix(0,reps,3)
tpI=matrix(0,reps,3)
for (rep in 1:reps){
  power[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  power[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,c(2,3)]=c(sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0),sum((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0))
  tpI[rep,1]=sum(((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),1])!=0)|((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[-c(2,3,4,5,6,7,10,11,12,13,14,15),2])!=0)) 
  print(c(rep,power[rep,]))
  print(c(rep,tpI[rep,]))
}
colSums(power)/(12*reps)
colSums(tpI)/(6*reps)

read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)

biassum=c(0,0,0)
rmsesum=c(0,0,0)
params=read.csv("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Para4.csv",row.names = 1)
Amat1=params[,1:2]
Dmat1=matrix(params[,7],20,(m-1))
para=cbind(Amat1,Dmat1)

for (rep in 1:reps){
  biassum[1:2]=biassum[1:2]+colSums(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)/10
  biassum[3]=biassum[3]+colMeans(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)
}
biassum/50

for (rep in 1:reps){
  rmsesum[1:2]=rmsesum[1:2]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,1:2]-Amat1)^2)
  rmsesum[3]=rmsesum[3]+colSums((read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,3]-Dmat1)^2)
}
sqrt(rmsesum[1:2]/(10*50))
sqrt(rmsesum[3]/(20*50))

absbias1=0
absbias2=0
for (rep in 1:reps){
  absbias1=absbias1+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),1]+0.5))
  absbias2=absbias2+sum(abs(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[c(2,3,4,5,6,7,10,11,12,13,14,15),2]+1))
}

absbias1/(12*50)
absbias2/(12*50)

###############
# Refit model #
###############
#grgamma=array(0,dim=c(r,r,J))
#grgamma[c(1,2),1,3:11]=matrix(t(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[1:9,]),2,9)
#grgamma[c(1,2),2,12:20]=matrix(t(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1)[10:18,]),2,9)
grgamma=matrix(0,J,y-1)
grgamma[3:20,]=as.matrix(read.csv(paste("/Users/zhux0445/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1))
grgamma[3:20,]=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1))

sparsity1=grgamma
for (j in 1:J){
  for (rr in 1:2){
    sparsity1[j,rr]=ifelse(grgamma[j,rr]==0,0,1)
  }
} 

sparsity2=grbeta
for (j in 1:J){
  for (rr in 1:2){
    sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
  }
} 

# use mirt for re-estimate
model3 <- '
    F1 = 1,3-11
    F2 = 2,12-20
    COV = F1*F2'
Group=c(rep('G1', N1), rep('G2', N2),rep('G3', N3))
#values <- multipleGroup(resp, model3, group = Group, pars = 'values')
#values
N1=N2=N3=1000 
N=N1+N2+N3
mirtCluster(8)

bias=matrix(0,reps,2)
rmse=matrix(0,reps,2)
absbias=matrix(0,reps,2)
for (rep in 1:reps){
  grgamma=matrix(0,J,y-1)
  grgamma[3:20,]=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1))
  grbeta=matrix(0,J,2)
  
  sparsity1=grgamma
  for (j in 1:J){
    for (rr in 1:2){
      sparsity1[j,rr]=ifelse(grgamma[j,rr]==0,0,1)
    }
  } 
  
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  } 

  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  Agp1=c(1,7,seq(11,51,5),seq(57,97,5))
  Agp23=rbind(c(106,211), c(112,217),c(116,221),c(121,226),c(126,231),c(131,236),c(136,241),c(141,246),c(146,251),c(151,256),c(156,261),c(162,267),c(167,272),c(172,277),c(177,282),c(182,287),c(187,292),c(192,297),c(197,302),c(202,307))
  Dgp1=seq(3,98,5)
  Dgp2=seq(108,203,5)
  Dgp3=seq(213,308,5)
  Dgp23=cbind(Dgp2,Dgp3)
  invsparsity1=1-sparsity1
  nonzerogamma=which(rowSums(invsparsity1)!=0)
  constrainA23<- cbind( Agp1,invsparsity1* Agp23)
  constrainA<- vector("list", length(nonzerogamma))
  for(nn in 1:length(nonzerogamma)){
    constrainA[[nn]]=c(fun.zero.omit((constrainA23[nonzerogamma,])[nn,]))
  }
  invsparsity2=1-sparsity2
  nonzerobeta=which(rowSums(invsparsity2)!=0)
  constrainD23<- cbind( Dgp1,invsparsity2* Dgp23)
  constrainD<- vector("list", length(nonzerobeta))
  for(nn in 1:length(nonzerobeta)){
    constrainD[[nn]]=c(fun.zero.omit((constrainD23[nonzerobeta,])[nn,]))
  }
  sim8reest <- multipleGroup(resp, model3, group = Group, constrain = c(constrainA,constrainD),invariance=c('free_means', 'free_var'))
 
  lowDIF=fun.zero.omit(c(coef(sim8reest,simplify=T)$G1$items[,1]-coef(sim8reest,simplify=T)$G2$items[,1],(coef(sim8reest,simplify=T)$G1$items[,2]-coef(sim8reest,simplify=T)$G2$items[,2])))
  highDIF=fun.zero.omit(c(coef(sim8reest,simplify=T)$G1$items[,1]-coef(sim8reest,simplify=T)$G3$items[,1],(coef(sim8reest,simplify=T)$G1$items[,2]-coef(sim8reest,simplify=T)$G3$items[,2])))
  bias[rep,1]=mean(lowDIF-0.5)
  bias[rep,2]=mean(highDIF-1)
  rmse[rep,1]=sqrt(sum((lowDIF-0.5)^2))
  rmse[rep,2]=sqrt(sum((highDIF-1)^2))
  absbias[rep,1]=mean(abs(lowDIF-0.5))
  absbias[rep,2]=mean(abs(highDIF-1))
  }


mean(na.omit(absbias[,1]))
mean(na.omit(absbias[,2]))


biasG=matrix(0,reps,1)
rmseG=matrix(0,reps,1)
absbiasG=matrix(0,reps,1)
biasAD=matrix(0,reps,3)
rmseAD=matrix(0,reps,3)

for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  StartVals=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1))
  grgamma=array(0,dim=c(y-1,r,J))
  grgamma[1,1,3:11]=t(as.matrix(StartVals[1:9,1]))
  grgamma[1,2,12:20]=t(as.matrix(StartVals[10:18,1]))
  grgamma
  grbeta=matrix(0,J,y-1)
  

gra=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,1:2])
grd=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,3])

mu100=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[1,]
mu200=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[2,]
#mu300=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[3,]
Sig100=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[4:5,]
Sig200=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[6:7,]
#Sig300=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[8:9,]
Mu.est=c(mu100,mu200)
Sig.est=rbind(Sig100,Sig200)
resp=resp[1:2000,]
if (min(resp)==0){
  resp2=as.matrix(resp)
  resp=resp+1
} else {
  resp2=as.matrix(resp)-1
}
N <- nrow(resp)
N.vec=c(1000,1000)
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

df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
iter <- 0
while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
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
  #LiA=E_step1( resp, N.vec, X, y, G, y.allgroup, Mu.list, Sig.list, gra, grd, grbeta, grgamma)
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
}
absbiasG[rep,]=mean(abs(fun.zero.omit(grgamma)+0.5))
biasG[rep,]=mean(fun.zero.omit(grgamma)+0.5)
rmseG[rep,]=sqrt(sum((fun.zero.omit(grgamma)+0.5)^2))
biasAD[rep,1:2]=colMeans(gra-Amat1)
rmseAD[rep,1:2]=sqrt(colSums((gra-Amat1)^2))
biasAD[rep,3]=colMeans(grd-Dmat1)
rmseAD[rep,3]=sqrt(colSums((grd-Dmat1)^2))
}

mean(na.omit(absbiasG))
colMeans(biasAD)
sqrt(colMeans(rmseAD^2))

for (rep in 1:reps){
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  StartVals=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Gamma8AdaptBIClowcor_",rep),row.names = 1))
  grgamma=array(0,dim=c(y-1,r,J))
  grgamma[1,1,3:11]=t(as.matrix(StartVals[1:9,2]))
  grgamma[1,2,12:20]=t(as.matrix(StartVals[10:18,2]))
  grbeta=matrix(0,J,y-1)
  
  
  gra=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,1:2])
  grd=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/ADmat8AdaptBIClowcor_",rep),row.names = 1)[,3])
  
  mu100=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[1,]
  mu200=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[2,]
  mu300=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[3,]
  Sig100=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[4:5,]
  Sig200=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[6:7,]
  Sig300=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/theta8AdaptBIClowcor_",rep),row.names = 1))[8:9,]
  Mu.est=c(mu100,mu300)
  Sig.est=rbind(Sig100,Sig300)
  resp=resp[c(1:1000,2001:3000),]
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  N.vec=c(1000,1000)
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
  
  df.a <- df.d  <- df.gamma <- df.beta <- df.Sig <- 1
  iter <- 0
  while(max(df.a)>eps | max(df.d)>eps | max(df.beta)>eps| max(df.gamma)>eps) # max(df.Mu)>eps | max(df.Sig)>eps |
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
    #LiA=E_step1( resp, N.vec, X, y, G, y.allgroup, Mu.list, Sig.list, gra, grd, grbeta, grgamma)
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
  }
  mean(abs(fun.zero.omit(grgamma)+1))
  gra-Amat
}

reps=38

bias2=matrix(0,reps,2)
rmse2=matrix(0,reps,2)
absbias2=matrix(0,reps,2)
for (rep in c(1:22,26:reps)){
  grgamma=matrix(0,J,y-1)
  grbeta=as.matrix(read.csv(paste("/Users/ruoyizhu/Documents/Github/RegDIF_SimData/Beta4adaptLowCor_",rep),row.names = 1))
  
  sparsity1=grgamma
  for (j in 1:J){
    for (rr in 1:2){
      sparsity1[j,rr]=ifelse(grgamma[j,rr]==0,0,1)
    }
  } 
  
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:2){
      sparsity2[j,rr]=ifelse(grbeta[j,rr]==0,0,1)
    }
  } 
  
  resp=responses[((rep-1)*N+1):((rep-1)*N+N1+N2+N3),]
  Agp1=c(1,7,seq(11,51,5),seq(57,97,5))
  Agp23=rbind(c(106,211), c(112,217),c(116,221),c(121,226),c(126,231),c(131,236),c(136,241),c(141,246),c(146,251),c(151,256),c(156,261),c(162,267),c(167,272),c(172,277),c(177,282),c(182,287),c(187,292),c(192,297),c(197,302),c(202,307))
  Dgp1=seq(3,98,5)
  Dgp2=seq(108,203,5)
  Dgp3=seq(213,308,5)
  Dgp23=cbind(Dgp2,Dgp3)
  invsparsity1=1-sparsity1
  nonzerogamma=which(rowSums(invsparsity1)!=0)
  constrainA23<- cbind( Agp1,invsparsity1* Agp23)
  constrainA<- vector("list", length(nonzerogamma))
  for(nn in 1:length(nonzerogamma)){
    constrainA[[nn]]=c(fun.zero.omit((constrainA23[nonzerogamma,])[nn,]))
  }
  invsparsity2=1-sparsity2
  nonzerobeta=which(rowSums(invsparsity2)!=0)
  constrainD23<- cbind( Dgp1,invsparsity2* Dgp23)
  constrainD<- vector("list", length(nonzerobeta))
  for(nn in 1:length(nonzerobeta)){
    constrainD[[nn]]=c(fun.zero.omit((constrainD23[nonzerobeta,])[nn,]))
  }
  sim4reest <- multipleGroup(resp, model3, group = Group, constrain = c(constrainA,constrainD),invariance=c('free_means', 'free_var'))
  
  lowDIF=fun.zero.omit(coef(sim4reest,simplify=T)$G1$items[,3]-coef(sim4reest,simplify=T)$G2$items[,3])
  highDIF=fun.zero.omit(coef(sim4reest,simplify=T)$G1$items[,3]-coef(sim4reest,simplify=T)$G3$items[,3])
  bias2[rep,1]=mean(lowDIF+0.5)
  bias2[rep,2]=mean(highDIF+1)
  rmse2[rep,1]=sqrt(sum((lowDIF+0.5)^2))
  rmse2[rep,2]=sqrt(sum((highDIF+1)^2))
  absbias2[rep,1]=mean(abs(lowDIF+0.5))
  absbias2[rep,2]=mean(abs(highDIF+1))
}


mean(na.omit(absbias2[,1]))
mean(na.omit(absbias2[,2]))

