sourceCpp("/Users/zhux0445/Documents/GitHub/RegDIF/Reg_DIF_EM.cpp")
sourceCpp("/Users/zhux0445/Documents/GitHub/RegDIF/matrix.cpp")
soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }


#E_step1=function(X=X,y=y,Mu.gp1=Mu.gp1,Mu.gp2=Mu.gp2,Mu.gp3=Mu.gp3,Sig.gp1=Sig.gp1,Sig.gp2=Sig.gp2,Sig.gp3=Sig.gp3,gra=gra,grgamma=grgamma){
E_step1=function( N.vec, X, y, y.allgroup, Mu.list, Sig.list, gra, grd, grbeta, grgamma){
  A.allgroups=numeric(nrow(X)*y) #A1, A2, A3
  for (yy in 1:y){
    A.allgroups[((yy-1)*nrow(X)+1):((yy-1)*nrow(X)+nrow(X))]=dmvnorm(X,Mu.list[((yy-1)*r+1):((yy-1)*r+r)],Sig.list[((yy-1)*r+1):((yy-1)*r+r),])
  }
  #calculation of n_g 
  axmat=gra%*%t(X) #a%*%X
  ygam.allgroups=matrix(0,J*y,(y-1)) #ygam1,ygam2,ygam3
  for (yy in 1:y){
    for (j in 1:J){
      ygam.allgroups[((yy-1)*J+j),]=y.allgroup[yy,]%*%grgamma[,,j]
    }
  }
  ygamat.allgrp=matrix(0,(J*y),nrow(X)) #ygamat1,ygamat2,ygamat3
  for (yy in 1:y){
    ygamat.allgrp[((yy-1)*J+1):((yy-1)*J+J),]=ygam.allgroups[((yy-1)*J+1):((yy-1)*J+J),]%*%t(X) 
  }
  grbeta.allgrp= matrix(0,(J*y),1) #grbeta1,grbeta2,grbeta3 "ncol=1" is the number of parameter beta, if we use multiple covariates later, change this value
  for (yy in 1:y){
    grbeta.allgrp[((yy-1)*J+1):((yy-1)*J+J),]=t(y.allgroup[yy,]%*%t(grbeta))
  }
  pstar.allgrp=array(double(J*(m-1)*G*y),dim = c(J,(m-1),G,y)) # pstar1,pstar2,pstar3
  p.allgrp=array(double(J*m*G*y),dim = c(J,m,G,y)) # p1,p2,p3
  for (yy in 1:y){
    for (g in 1:G)
    {
      pstar.allgrp[,,g,yy] = 1/(1+exp(-(grd+axmat[,g]+ygamat.allgrp[((yy-1)*J+1):((yy-1)*J+J),g]+grbeta.allgrp[((yy-1)*J+1):((yy-1)*J+J),])))
      p.allgrp[,,g,yy] = t(-diff(rbind(rep(1,J),t(pstar.allgrp[,,g,yy]),rep(0,J))))
    }
  }
  
  LiA=matrix(double(N*G),N,G)
  for (yy in 1:y){
    pij=array(double(J*G*N.vec[yy]),dim=c(J,G,N.vec[yy]))
    for (j in 1:J)
    {
      for (g in 1:G)
      {
        pij[j,g,]=p.allgrp[j,(resp[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy])),j]),g,yy]
      }
    }
    LiA[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy])),]=t(apply(pij, c(2,3), prod)* A.allgroups[((yy-1)*nrow(X)+1):((yy-1)*nrow(X)+nrow(X))])
  }
  return(LiA)
}
  
ng.est= function(LiA,y,N.vec){
  Pi = apply(LiA,1,sum)
  ng = apply(LiA/Pi,2,sum)
  ng.allgrp=numeric(G*y)
  for (yy in 1:y){
    ng.allgrp[((yy-1)*G+1):((yy-1)*G+G)]=apply(LiA[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy])),]/Pi[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy]))],2,sum)
  }
  return(c(ng,ng.allgrp))
}

rgk.est=function(j,Xijk,LiA,y,N.vec){
  Pi = apply(LiA,1,sum)
  rLiA <- array(double(N*G*m),dim = c(N,G,m))
  for(k in 1:m){
    rLiA[,,k]=Xijk[,j,k]*LiA/Pi
  }
  rgk <- apply(rLiA,c(2,3),sum)
  rgk.allgrp=matrix(0,G*y,m)
  for (yy in 1:y){
    rgk.allgrp[((yy-1)*G+1):((yy-1)*G+G),]= apply(rLiA[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy])),,],c(2,3),sum)
  }
  return(rbind(rgk,rgk.allgrp))
}

M_step=function(j,grd,gra,grgamma,grbeta,max.tol,X,y.allgroup,y,G,eta){
  d <- grd[j,] 
  a <- gra[j,]
  gam=grgamma[,,j]
  bet=grbeta[j,]
  #M-step loop starts for item j
  miter <- 0
  add <- max.tol+1
  while(sum(abs(add))>max.tol)
  {
    miter <- miter+1
    for  (yy in 1:y){
      for(g in 1:G){
        Pstar[g,,yy] <- 1/(1+exp(-(d+rep(a%*%X[g,])+y.allgroup[yy,]%*%bet+rep((y.allgroup[yy,]%*%gam)%*%X[g,]))))
        P[g,,yy] <- -diff(c(1,Pstar[g,,yy],0))
      }
    }
    Qstar <- 1-Pstar
    
    #calculating the score vector
    if (m==2){
      Dsco=numeric(m-1)
      for (yy in 1:y){
        Dsco <- Dsco+sum(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy],1,diff)*Pstar[,,yy]*Qstar[,,yy])
      }
    } else {
      Dsco=numeric(m-1)
      for (yy in 1:y){
        Dsco <- Dsco+apply(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy],1,diff)*Pstar[,,yy]*Qstar[,,yy],2,sum)
      }
    }
    PQdif=array(double(G*m*y),dim=c(G,m,y))
    for (yy in 1:y){
      PQdif[,,yy] <- -t(apply(cbind(0,Pstar[,,yy]*Qstar[,,yy],0),1, diff))
    }
    Asco = numeric(sum(ifelse(a==0,0,1)))
    for (yy in 1:y){
      Asco=Asco+(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy]*PQdif[,,yy],1,sum))%*%(X%*%ifelse(a==0,0,1))
    }
    Gamsco =numeric(sum(ifelse(gam==0,0,1)))
    if (sum(ifelse(gam==0,0,1))>0){
      Gamsco.all=gam
      for (yy in 2:y){
        Gamsco.all[yy-1,]=(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy]*PQdif[,,yy],1,sum))%*%X
      }
      Gamsco=as.vector(Gamsco.all)[which(gam!=0)]
    } else {
      Gamsco=as.null(Gamsco)
    }
    Betsco =numeric(sum(ifelse(bet==0,0,1)))
    if (sum(ifelse(bet==0,0,1))>0){
      Betsco.all=bet
      if (m==2){
        for (yy in 2:y){
          Betsco.all[yy-1] <- sum(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy],1,diff)*Pstar[,,yy]*Qstar[,,yy])
        }
        Betsco <- Betsco.all[which(bet!=0)]
      } else {
        for (yy in 2:y){
          Betsco.all[yy-1] <- apply(apply(rgk[(yy*G+1):(yy*G+G),]/P[,,yy],1,diff)*Pstar[,,yy]*Qstar[,,yy],2,sum)
        }
        Betsco <- Betsco.all[which(bet!=0)]
      }
        
    } else {
      Betsco=as.null(Betsco)
    }
    
      minusgrad <- -c(Dsco,Asco, Gamsco, Betsco)
      grad <- c(Dsco,Asco,Gamsco, Betsco)
      FI <- matrix(0,length(minusgrad),length(minusgrad))
      for (kk in 1:length(Dsco)){
        for (yy in 1:y){
          FI[kk,kk] = FI[kk,kk]+ (-sum(ng[(yy*G+1):(yy*G+G)]*(Pstar[,kk,yy]*Qstar[,kk,yy])^2*(1/P[,kk,yy]+1/P[,kk+1,yy])))
        }
      }
      if (m>2){
        for (mm in 2:(m-1)){
          FI[mm,mm-1] <-FI[mm-1,mm] <- sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm-1]*Qstar[,mm-1]*(1/P[,mm]))#(2,1),(3,2)
        }
        #for (mm in 1:(m-2)){
        #  FI[mm,mm+1] <-  sum(ng*Pstar[,mm]*Qstar[,mm]*Pstar[,mm+1]*Qstar[,mm+1]*(1/P[,mm+1]))#(1,2),(2,3)
        #}
      }
      for (kk in (length(Dsco)+1):(length(Dsco)+length(Asco))){
        for (yy in 1:y){
          FI[kk,kk] = FI[kk,kk]+ (-sum(ng[(yy*G+1):(yy*G+G)]*(X%*%ifelse(a==0,0,1))[,(kk-length(Dsco))]*(X%*%ifelse(a==0,0,1))[,(kk-length(Dsco))]*apply(PQdif[,,yy]^2/P[,,yy],1,sum)))
          for (bb in 1:length(Dsco)){
            FI[kk,bb] <-  FI[kk,bb]+ sum(ng[(yy*G+1):(yy*G+G)]*(X%*%ifelse(a==0,0,1))[,(kk-length(Dsco))]*(Pstar[,bb,yy]*Qstar[,bb,yy])*(PQdif[,bb,yy]/P[,bb,yy]-PQdif[,bb+1,yy]/P[,bb+1,yy]))
            FI[bb,kk] <-  FI[kk,bb]
          }
        }
       
        
        #if (length(Asco)>1){
        #  n.cross=choose(length(Asco),2)
        #  for (c in 1:n.cross){
        #    for (yy in 1:y){
        #      FI[kk,kk] = FI[kk,kk]+ (-sum(ng[(yy*G+1):(yy*G+G)]*(X%*%ifelse(a==0,0,1))[,(kk-length(Dsco))]*(X%*%ifelse(a==0,0,1))[,(kk-length(Dsco))]*apply(PQdif[,,yy]^2/P[,,yy],1,sum)))
        #    }
        #  }
        #}
      }
     
      if (length(Gamsco)>0){
        for (kk in (length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))){
          grp.number=which(gam!=0)[kk-(length(Dsco)+length(Asco))]%%(y-1)
          if (grp.number==0){
            grp=y
            dim.number=which(gam!=0)[kk-(length(Dsco)+length(Asco))]/(y-1)
          } else {
            grp=grp.number+1
            dim.number=which(gam!=0)[kk-(length(Dsco)+length(Asco))]%/%(y-1)+1
          }
          FI[kk,kk] = FI[kk,2]= FI[2,kk]= -sum(ng[(grp*G+1):(grp*G+G)]*X[,dim.number]*X[,dim.number]*apply(PQdif[,,grp]^2/P[,,grp],1,sum))
          #if (length(Asco)>1){
          #  for (aa in 1:length(asco)){
          #    FI[kk,aa] = FI[aa,kk]=-sum(ng[(grp*G+1):(grp*G+G)]*X[,dim.number]*X[,aa]*apply(PQdif[,,grp]^2/P[,,grp],1,sum))
          #  }
          #}
          for (dd in 1:length(Dsco)){
            FI[kk,dd]=FI[dd,kk]=sum(ng[(grp*G+1):(grp*G+G)]*X[,dim.number]*Pstar[,dd,grp]*Qstar[,dd,grp]*(PQdif[,dd,grp]/P[,dd,grp]-PQdif[,dd+1,grp]/P[,dd+1,grp]))
          }
        }
      }
    
    
    add <- qr.solve(FI,minusgrad)
    d=d+add[1:(m-1)]
    a[which(a!=0)]= a[which(a!=0)]+add[(length(Dsco)+1):(length(Dsco)+length(Asco))]
    if (length(Gamsco)>0){
      gam0=gam
      gam0[which(gam!=0)]=gam[which(gam!=0)]+add[(length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))]
      for (mm in (length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))){
        add[mm]=soft(  (gam0[which(gam!=0)])[mm-(length(Dsco)+length(Asco))],-eta/FI[mm,mm])- (gam[which(gam!=0)])[mm-(length(Dsco)+length(Asco))]
      }
      gam[which(gam!=0)]=gam[which(gam!=0)]+add[(length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))]
      
    }
  }
  return(c(d=d,a=a,gam=gam,bet=bet))
  #end of M step loop
}

sumoverk=function(G,rgky,aj,dj,gamjy,X){
  sumoverky=numeric(G)
  for(g in 1:G){
    sumoverky[g]=rgky[g,]%*%log(-diff(c(1,1/(1+exp(-(dj+rep(aj%*%X[g,])+rep(gamjy%*%X[g,])))),0)))
  }
  return(sum(sumoverky))
}
?dmvnorm
