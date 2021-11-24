#MBR three DIF functions
library(MASS)
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)

sourceCpp("Reg_DIF_EM.cpp")
soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

M_step=function(j,ng,rgk,grd,gra,grgamma,grbeta,max.tol,X,y.allgroup,y,G,m,eta){
  d <- grd[j,] 
  a <- gra[j,]
  gam=grgamma[,,j]
  bet=grbeta[j,]
  Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
  P<- array(double(G*m*(y)),dim=c(G,m,y))
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
      Gamsco.all=matrix(gam,y-1,r)
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
          FI[kk,dd]=FI[dd,kk]=sum(ng[(grp*G+1):(grp*G+G)]*X[,dim.number]*Pstar[,dd,grp]*Qstar[,dd,grp]*(PQdif[,dd,grp]/P[,dd,grp]-PQdif[,dd+1,grp]/P[,dd+1,grp])) #GRM
        }
      }
    }
    
    if (length(Betsco)>0){
      for (kk in (length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))){
        grp.number=which(bet!=0)[kk-(length(Dsco)+length(Asco)+length(Gamsco))]
        FI[kk,kk] =  FI[kk,m-1] =  FI[m-1,kk]= (-sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1])^2*(1/P[,1,grp.number+1]+1/P[,m,grp.number+1])))
        FI[kk,m] =  FI[m,kk]= sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(X%*%ifelse(a==0,0,1))*Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1]*(PQdif[,1,grp.number+1]/P[,1,grp.number+1]-PQdif[,2,grp.number+1]/P[,2,grp.number+1]))
        if (length(Gamsco)>0){
          for (kk2 in (length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))){
            grp.number2=which(gam!=0)[kk2-(length(Dsco)+length(Asco))]%%(y-1)
            if (grp.number2==0){
              grp=y
              dim.number=which(gam!=0)[kk2-(length(Dsco)+length(Asco))]/(y-1)
            } else {
              grp=grp.number2+1
              dim.number=which(gam!=0)[kk2-(length(Dsco)+length(Asco))]%/%(y-1)+1
            }
            if (grp==(grp.number+1)){
              FI[kk,kk2] =  FI[kk2,kk]= sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(X%*%ifelse(a==0,0,1))*Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1]*(PQdif[,1,grp.number+1]/P[,1,grp.number+1]-PQdif[,2,grp.number+1]/P[,2,grp.number+1])) #2pl only, not for GRM
            }
          }
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
    if (length(Betsco)>0){
      bet0=bet
      bet0[which(bet!=0)]=bet[which(bet!=0)]+add[(length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))]
      for (mm in (length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))){
        add[mm]=soft(  (bet0[which(bet!=0)])[mm-(length(Dsco)+length(Asco)+length(Gamsco))],-eta/FI[mm,mm])- (bet[which(bet!=0)])[mm-(length(Dsco)+length(Asco)+length(Gamsco))]
      }
      bet[which(bet!=0)]=bet[which(bet!=0)]+add[(length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))]
      
    }
    #print(c(miter,add))
  }
  return(c(d=d,a=a,gam=gam,bet=bet))
  #end of M step loop
}

M_step_Adaptive=function(j,ng,rgk,grd,gra,grgamma,grgamma00,grbeta,grbeta00,max.tol,X,y.allgroup,y,G,m,eta,lam){
  d <- grd[j,] 
  a <- gra[j,]
  gam=grgamma[,,j]
  gammle=grgamma00[,,j]
  bet=grbeta[j,]
  betmle=grbeta00[j,]
  Pstar <- Qstar <- array(double(G*(m-1)*(y)),dim=c(G,m-1,y))
  P<- array(double(G*m*(y)),dim=c(G,m,y))
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
    
    if (length(Betsco)>0){
      for (kk in (length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))){
        grp.number=which(bet!=0)[kk-(length(Dsco)+length(Asco)+length(Gamsco))]
        FI[kk,kk] =  FI[kk,m-1] =  FI[m-1,kk]= (-sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1])^2*(1/P[,1,grp.number+1]+1/P[,m,grp.number+1])))
        FI[kk,m] =  FI[m,kk]= sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(X%*%ifelse(a==0,0,1))*Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1]*(PQdif[,1,grp.number+1]/P[,1,grp.number+1]-PQdif[,2,grp.number+1]/P[,2,grp.number+1]))
        if (length(Gamsco)>0){
          for (kk2 in (length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))){
            grp.number2=which(gam!=0)[kk2-(length(Dsco)+length(Asco))]%%(y-1)
            if (grp.number2==0){
              grp=y
              dim.number=which(gam!=0)[kk2-(length(Dsco)+length(Asco))]/(y-1)
            } else {
              grp=grp.number2+1
              dim.number=which(gam!=0)[kk2-(length(Dsco)+length(Asco))]%/%(y-1)+1
            }
            if (grp==(grp.number+1)){
              FI[kk,kk2] =  FI[kk2,kk]= sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(X%*%ifelse(a==0,0,1))*Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1]*(PQdif[,1,grp.number+1]/P[,1,grp.number+1]-PQdif[,2,grp.number+1]/P[,2,grp.number+1])) #2pl only, not for GRM
            }
          }
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
        add[mm]=soft(  (gam0[which(gam!=0)])[mm-(length(Dsco)+length(Asco))],-(eta/(abs((gammle[which(gam!=0)])[mm-(length(Dsco)+length(Asco))])^(lam)))/FI[mm,mm])- (gam[which(gam!=0)])[mm-(length(Dsco)+length(Asco))]
      }
      gam[which(gam!=0)]=gam[which(gam!=0)]+add[(length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))]
      
    }
    if (length(Betsco)>0){
      bet0=bet
      bet0[which(bet!=0)]=bet[which(bet!=0)]+add[(length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))]
      for (mm in (length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))){
        add[mm]=soft(  (bet0[which(bet!=0)])[mm-(length(Dsco)+length(Asco)+length(Gamsco))],-(eta/(abs((betmle[which(bet!=0)])[mm-(length(Dsco)+length(Asco)+length(Gamsco))])^(lam)))/FI[mm,mm])- (bet[which(bet!=0)])[mm-(length(Dsco)+length(Asco)+length(Gamsco))]
      }
      bet[which(bet!=0)]=bet[which(bet!=0)]+add[(length(Dsco)+length(Asco)+length(Gamsco)+1):(length(Dsco)+length(Asco)+length(Gamsco)+length(Betsco))]
      
    }
  }
  return(c(d=d,a=a,gam=gam,bet=bet))
  #end of M step loop
}


#Starting Values

  DIF_init=function(u,Group,indic){
    m=2 ##fixed, 2pl only
    N=nrow(u)
    J=ncol(u)
    domain=nrow(indic)
    y=length(unique(Group)) 
    y.allgroup=rbind(rep(0,y-1),diag(y-1)) 
    G=matrix(0,N,y-1)
    for (yy in 1:y){
      vec=which(Group==sort(unique(Group))[yy])
      for (i in 1:length(vec)){
        G[vec[i],]=y.allgroup[yy,]
      }
    }
    COV <- matrix(TRUE,domain,domain); diag(COV)=FALSE
    model <- mirt.model(t(indic), COV=COV) ##
    md.noncons0 <- multipleGroup(u, model, group = Group,SE=TRUE)
    starting5new=cbind(coef(md.noncons0,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,1:domain]-coef(md.noncons0,simplify=T)[[1]]$items[,1:domain])
    }
    for (yy in 2:y){
      starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,domain+1]-coef(md.noncons0,simplify=T)[[1]]$items[,domain+1])
      
    }
    gra00=as.matrix(starting5new[,1:domain])
    rownames(gra00) <- c()
    colnames(gra00) <- c()
    grd00=matrix(starting5new[,domain+1],J,1)
    grgamma00=array(0,dim=c((y-1),domain,J))
    for (yy in 2:y){
      grgamma00[(yy-1),,]=t(starting5new[,(domain+2+(yy-2)*domain):(2*domain+1+(yy-2)*domain)])
    }
    grbeta00=t(as.matrix(starting5new[,(2*domain+2+(y-2)*domain):(2*domain+1+(y-2)*domain+m)]))
    rownames(grbeta00) <- c()
    colnames(grbeta00) <- c()
    Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
    Mu0 = matrix(0,domain,y)
    for (yy in 1:y){
      Sigma0[,,yy]=coef(md.noncons0,simplify=T)[[yy]]$cov
      #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
    }
    return(list(G=G,y=y,r=domain,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Sigma0=Sigma0,Mu0 =Mu0))
  }

##################################################
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

Reg_DIF <- function(resp,Group,indic,eta,eps =1e-3,max.tol=1e-7)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  #m,r,y,N.vec,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL
  init=DIF_init(u=resp2,Group=Group,indic=indic)
  m=2 #fixed 2pl
  r=init$r
  y=init$y
  N.vec=as.vector(table(Group))
  gra00=init$gra00
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
  grgamma=array(0,dim=c((y-1),r,J));grgamma[1,1,4:5]=grgamma[1,2,12:13]=-0.5;grgamma[2,1,4:5]=grgamma[2,2,12:13]=-1
  #grgamma=grgamma00 #array(0,dim=c((y-1),r,J))
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


#####shiny app######
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for data upload app ----
ui <- navbarPage("RegDIF Shiny App",
                 tabPanel("Lasso EM",
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              ########################################################
                              # Input: Select a file ----
                              fileInput("file1", "Choose Data CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header1", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep1", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # # Input: Select number of rows to display ----
                              radioButtons("disp1", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: number of response categories
                              tags$hr(),
                              numericInput("m", "number of response categories", min=2, max=10, value=2),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: domain
                              tags$hr(),
                              numericInput("r", "domain", min=1, max=3, value=2),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: number of groups
                              tags$hr(),
                              numericInput("y", "number of groups", min=2, max=10, value=3),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: N.vec
                              fileInput("file2", "Choose Group Sample Size CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              helpText("Note: The Group Sample Size file should be a vector including sample size of each covariate group."),
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header2", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep2", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              ########################################################
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: eta
                              tags$hr(),
                              numericInput("eta", "Tuning Parameter", min=0, max=100, value=0),
                              
                              ########################################################
                              #Horizontal line ----
                              tags$hr(),
                              
                              # Input: Starting Values for parameter a
                              fileInput("file3", "Choose Starting Values for parameter a CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header3", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep3", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp3", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: Starting Values for parameter d
                              fileInput("file4", "Choose Starting Values for parameter d CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header4", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep4", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp4", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: Starting Values for parameter gamma
                              fileInput("file5", "Choose Starting Values for DIF on slope CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header5", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep5", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp5", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              #Horizontal line ----
                              tags$hr(),
                              ########################################################
                              # Input: Starting Values for parameter beta
                              fileInput("file6", "Choose Starting Values for DIF on intercept CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header6", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep6", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp6", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              # Horizontal line ----
                              tags$hr(),
                              
                              ########################################################
                              # Input: Starting Values for parameter mu
                              fileInput("file7", "Choose Starting Values for mean vectors CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header7", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep7", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp7", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              # Horizontal line ----
                              tags$hr(),
                              
                              ########################################################
                              # Input: Starting Values for parameter sigma
                              fileInput("file8", "Choose Starting Values for covariance matrix CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Checkbox if file has header ----
                              checkboxInput("header8", "Header", TRUE),
                              
                              # Input: Select separator ----
                              radioButtons("sep8", "Separator",
                                           choices = c(Comma = ",",
                                                       Semicolon = ";",
                                                       Tab = "\t"),
                                           selected = ","),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              #Input: Select number of rows to display ----
                              radioButtons("disp8", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go3", "Run"),
                              
                              tags$hr(),
                              radioButtons("checkGroup2", "Download Results",
                                           choices = list("All Results" = "all",
                                                          "Item Parameters" = "item",
                                                          "Covariance Matrix" = "cov"),
                                           selected = "all"),
                              #Button
                              downloadButton("downloadDataall2", "Download Results")
                              #downloadButton("downloadDataall2", "Download All Results")
                              # downloadButton("downloadDataitem2", "Download Item Parameter Results"),
                              # downloadButton("downloadDatasigma2", "Download Covariance Results")
                              
                              
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents1"),
                              
                              
                              h2("Group Sample Size Matrix"),
                              #Output: N.vec,
                              tableOutput("contents2"),
                              
                              h2("Starting Values for parameter a"),
                              #Output: gra00,
                              tableOutput("contents3"),
                              
                              h2("Starting Values for parameter d"),
                              #Output: grd00,
                              tableOutput("contents4"),
                              
                              h2("Starting Values for DIF on slope"),
                              #Output: grgamma00,
                              tableOutput("contents5"),
                              
                              h2("Starting Values for DIF on intercept"),
                              #Output: grbeta00,
                              tableOutput("contents6"),
                              
                              h2("Starting Values for mean vector"),
                              #Output: grgamma00,
                              tableOutput("contents7"),
                              
                              h2("Starting Values for covariance matrix"),
                              #Output: grbeta00,
                              tableOutput("contents8"),
                              
                              
                              h2("Item Parameter Results"),
                              tableOutput("par1"),
                              
                              h2("Covariance Matrix"),
                              tableOutput("cov1"),
                              
                              h2("Model Fit"),
                              textOutput("bic1"),
                              
                              
                            )
                            
                          )
                 ),tabPanel("Lasso EMM",
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                ########################################################
                                # Input: Select a file ----
                                fileInput("file1", "Choose Data CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header1", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep1", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # # Input: Select number of rows to display ----
                                radioButtons("disp1", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: number of response categories
                                tags$hr(),
                                numericInput("m", "number of response categories", min=2, max=10, value=2),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: domain
                                tags$hr(),
                                numericInput("r", "domain", min=1, max=3, value=2),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: number of groups
                                tags$hr(),
                                numericInput("y", "number of groups", min=2, max=10, value=3),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: N.vec
                                fileInput("file2", "Choose Group Sample Size CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                helpText("Note: The Group Sample Size file should be a vector including sample size of each covariate group."),
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header2", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep2", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                ########################################################
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: eta
                                tags$hr(),
                                numericInput("eta", "Tuning Parameter", min=0, max=100, value=0),
                                
                                ########################################################
                                #Horizontal line ----
                                tags$hr(),
                                
                                # Input: Starting Values for parameter a
                                fileInput("file3", "Choose Starting Values for parameter a CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header3", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep3", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp3", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: Starting Values for parameter d
                                fileInput("file4", "Choose Starting Values for parameter d CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header4", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep4", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp4", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: Starting Values for parameter gamma
                                fileInput("file5", "Choose Starting Values for DIF on slope CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header5", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep5", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp5", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: Starting Values for parameter beta
                                fileInput("file6", "Choose Starting Values for DIF on intercept CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header6", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep6", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp6", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                # Horizontal line ----
                                tags$hr(),
                                
                                ########################################################
                                # Input: Starting Values for parameter mu
                                fileInput("file7", "Choose Starting Values for mean vectors CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header7", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep7", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp7", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                # Horizontal line ----
                                tags$hr(),
                                
                                ########################################################
                                # Input: Starting Values for parameter beta
                                fileInput("file8", "Choose Starting Values for covariance matrix CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header8", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep8", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp8", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                # Horizontal line ----
                                tags$hr(),
                                actionButton("go3", "Run"),
                                
                                tags$hr(),
                                radioButtons("checkGroup2", "Download Results",
                                             choices = list("All Results" = "all",
                                                            "Item Parameters" = "item",
                                                            "Covariance Matrix" = "cov"),
                                             selected = "all"),
                                #Button
                                downloadButton("downloadDataall2", "Download Results")
                                #downloadButton("downloadDataall2", "Download All Results")
                                # downloadButton("downloadDataitem2", "Download Item Parameter Results"),
                                # downloadButton("downloadDatasigma2", "Download Covariance Results")
                                
                                
                                
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(
                                
                                h2("Data"),
                                # Output: Data file ----
                                tableOutput("contents1"),
                                
                                
                                h2("Group Sample Size Matrix"),
                                #Output: N.vec,
                                tableOutput("contents2"),
                                
                                h2("Starting Values for parameter a"),
                                #Output: gra00,
                                tableOutput("contents3"),
                                
                                h2("Starting Values for parameter d"),
                                #Output: grd00,
                                tableOutput("contents4"),
                                
                                h2("Starting Values for DIF on slope"),
                                #Output: grgamma00,
                                tableOutput("contents5"),
                                
                                h2("Starting Values for DIF on intercept"),
                                #Output: grbeta00,
                                tableOutput("contents6"),
                                
                                h2("Starting Values for mean vector"),
                                #Output: grgamma00,
                                tableOutput("contents7"),
                                
                                h2("Starting Values for covariance matrix"),
                                #Output: grbeta00,
                                tableOutput("contents8"),
                                
                                
                                h2("Item Parameter Results"),
                                tableOutput("par1"),
                                
                                h2("Covariance Matrix"),
                                tableOutput("cov1"),
                                
                                h2("Model Fit"),
                                textOutput("bic1"),
                                
                                
                              )
                              
                            )
                 ),tabPanel("Adaptive Lasso",
                            # Sidebar layout with input and output definitions ----
                            sidebarLayout(
                              
                              # Sidebar panel for inputs ----
                              sidebarPanel(
                                ########################################################
                                # Input: Select a file ----
                                fileInput("file1", "Choose Data CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header1", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep1", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # # Input: Select number of rows to display ----
                                radioButtons("disp1", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: number of response categories
                                tags$hr(),
                                numericInput("m", "number of response categories", min=2, max=10, value=2),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: domain
                                tags$hr(),
                                numericInput("r", "domain", min=1, max=3, value=2),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: number of groups
                                tags$hr(),
                                numericInput("y", "number of groups", min=2, max=10, value=3),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: N.vec
                                fileInput("file2", "Choose Group Sample Size CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                helpText("Note: The Group Sample Size file should be a vector including sample size of each covariate group."),
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header2", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep2", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                ########################################################
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: eta
                                tags$hr(),
                                numericInput("eta", "Tuning Parameter", min=0, max=100, value=0),
                                
                                ########################################################
                                #Horizontal line ----
                                tags$hr(),
                                
                                # Input: Starting Values for parameter a
                                fileInput("file3", "Choose Starting Values for parameter a CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header3", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep3", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp3", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: Starting Values for parameter d
                                fileInput("file4", "Choose Starting Values for parameter d CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header4", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep4", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp4", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: Starting Values for parameter gamma
                                fileInput("file5", "Choose Starting Values for DIF on slope CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header5", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep5", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp5", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                
                                #Horizontal line ----
                                tags$hr(),
                                ########################################################
                                # Input: Starting Values for parameter beta
                                fileInput("file6", "Choose Starting Values for DIF on intercept CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header6", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep6", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp6", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                # Horizontal line ----
                                tags$hr(),
                                
                                ########################################################
                                # Input: Starting Values for parameter mu
                                fileInput("file7", "Choose Starting Values for mean vectors CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header7", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep7", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp7", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                # Horizontal line ----
                                tags$hr(),
                                
                                ########################################################
                                # Input: Starting Values for parameter beta
                                fileInput("file8", "Choose Starting Values for covariance matrix CSV File",
                                          multiple = TRUE,
                                          accept = c("text/csv",
                                                     "text/comma-separated-values,text/plain",
                                                     ".csv")),
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                # Input: Checkbox if file has header ----
                                checkboxInput("header8", "Header", TRUE),
                                
                                # Input: Select separator ----
                                radioButtons("sep8", "Separator",
                                             choices = c(Comma = ",",
                                                         Semicolon = ";",
                                                         Tab = "\t"),
                                             selected = ","),
                                
                                
                                # Horizontal line ----
                                tags$hr(),
                                
                                #Input: Select number of rows to display ----
                                radioButtons("disp8", "Display",
                                             choices = c(Head = "head",
                                                         All = "all"),
                                             selected = "head"),
                                # Horizontal line ----
                                tags$hr(),
                                actionButton("go3", "Run"),
                                
                                tags$hr(),
                                radioButtons("checkGroup2", "Download Results",
                                             choices = list("All Results" = "all",
                                                            "Item Parameters" = "item",
                                                            "Covariance Matrix" = "cov"),
                                             selected = "all"),
                                #Button
                                downloadButton("downloadDataall2", "Download Results")
                                #downloadButton("downloadDataall2", "Download All Results")
                                # downloadButton("downloadDataitem2", "Download Item Parameter Results"),
                                # downloadButton("downloadDatasigma2", "Download Covariance Results")
                                
                                
                                
                              ),
                              
                              # Main panel for displaying outputs ----
                              mainPanel(
                                
                                h2("Data"),
                                # Output: Data file ----
                                tableOutput("contents1"),
                                
                                
                                h2("Group Sample Size Matrix"),
                                #Output: N.vec,
                                tableOutput("contents2"),
                                
                                h2("Starting Values for parameter a"),
                                #Output: gra00,
                                tableOutput("contents3"),
                                
                                h2("Starting Values for parameter d"),
                                #Output: grd00,
                                tableOutput("contents4"),
                                
                                h2("Starting Values for DIF on slope"),
                                #Output: grgamma00,
                                tableOutput("contents5"),
                                
                                h2("Starting Values for DIF on intercept"),
                                #Output: grbeta00,
                                tableOutput("contents6"),
                                
                                h2("Starting Values for mean vector"),
                                #Output: grgamma00,
                                tableOutput("contents7"),
                                
                                h2("Starting Values for covariance matrix"),
                                #Output: grbeta00,
                                tableOutput("contents8"),
                                
                                
                                h2("Item Parameter Results"),
                                tableOutput("par1"),
                                
                                h2("Covariance Matrix"),
                                tableOutput("cov1"),
                                
                                h2("Model Fit"),
                                textOutput("bic1"),
                                
                                
                              )
                              
                            )
                 ))

# Define server logic to read selected file ----
server <- function(input, output,session) {
  url2 <- a("M2PL estimation", href="https://www.google.com/",)
  output$tab2 <- renderUI({
    tagList(url2)
  })
  
  #Lasso-EM output
  resp1<-reactive({req(input$file1)
    df <- data.matrix(read.csv(input$file1$datapath,
                               header = input$header1,
                               sep = input$sep1))
    return(df)
  })
  output$contents1 <- renderTable({
    if(input$disp1 == "head") {
      return(resp1()[1:6,])
    }
    else {
      return(resp1())
    }
    
  })
  ##############################################################
  N.vec1<-reactive({req(input$file2)
    df2 <- data.matrix(read.csv(input$file2$datapath,
                                header = input$header2,
                                sep = input$sep2))
    return(df2)
  })
  
  output$contents2 <- renderTable({

      return(N.vec1())

  })
  ##############################################################
  m1<-reactive({req(input$m)})
  
  r1<-reactive({req(input$r)})
  y1<-reactive({req(input$y)})
  eta1<-reactive({req(input$eta)})
  ##############################################################
  gra001<-reactive({req(input$file3)
    df3 <- data.matrix(read.csv(input$file3$datapath,
                               header = input$header3,
                               sep = input$sep3))
    return(df3)
  })
  output$contents3 <- renderTable({
    if(input$disp3 == "head") {
      return(gra001()[1:6,])
    }
    else {
      return(gra001())
    }
    
  })
  ##############################################################
  grd001<-reactive({req(input$file4)
    df4 <- data.matrix(read.csv(input$file4$datapath,
                                header = input$header4,
                                sep = input$sep4))
    return(df4)
  })
  output$contents4 <- renderTable({
    if(input$disp4 == "head") {
      return(grd001()[1:6,])
    }
    else {
      return(grd001())
    }
    
  })
  
  ##############################################################
  grgamma001<-reactive({req(input$file5)
    df5 <- data.matrix(read.csv(input$file5$datapath,
                                header = input$header5,
                                sep = input$sep5))
    return(df5)
  })
  output$contents5 <- renderTable({
    if(input$disp5 == "head") {
      return(grgamma001()[1:6,])
    }
    else {
      return(grgamma001())
    }
    
  })
  
  ##############################################################
  grbeta001<-reactive({req(input$file6)
    df6 <- data.matrix(read.csv(input$file6$datapath,
                                header = input$header6,
                                sep = input$sep6))
    return(df6)
  })
  output$contents6 <- renderTable({
    if(input$disp6 == "head") {
      return(grbeta001()[1:6,])
    }
    else {
      return(grbeta001())
    }
    
  })
  ##############################################################
  Mu.list1<-reactive({req(input$file7)
    df7 <- data.matrix(read.csv(input$file7$datapath,
                                header = input$header7,
                                sep = input$sep7))
    return(df7)
  })
  output$contents7 <- renderTable({
    return( Mu.list1())
  })
  ##############################################################
  Sig.list1<-reactive({req(input$file8)
    df8 <- data.matrix(read.csv(input$file8$datapath,
                                header = input$header8,
                                sep = input$sep8))
    return(df8)
  })
  output$contents8 <- renderTable({
    return(Sig.list1())
  })
  ##############################################################
  result<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    result<-Reg_DIF(resp1(),m1(),r1(),y1(),N.vec1(),eta1(),1e-3,1e-7,gra001(),grd001(),grbeta001(),grgamma001(),Mu.list1(),Sig.list1()) 
    return(result)
  })
  
  output$par2<-renderTable({
    input$go3
    isolate({
     # m<-cbind(result()$ra,result()$rb)
      m<-cbind(result()$Beta)
      #colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
      return(m)
    })
    
  })
  
  output$cov2<-renderTable({
    input$go3
    isolate({
      m<-result()$Covs
      #rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
      return(m)
    })
    
  })
  
  #output$gic2<-renderText({
   # input$go3
  #  isolate({
  #    paste("GIC =", round(result()$GIC,2))
  #  })
  #})
  
  #output$aic2<-renderText({
  #  input$go3
  #  isolate({
  #    paste("AIC =", round(result()$AIC,2))
  #  })
  #})
  
  output$bic2<-renderText({
    input$go3
    isolate({
      paste("BIC =", round(result()$bic,2))
    })
  })
  #Downloadable csv of selected dataset ----
  output$downloadDataall2 <- downloadHandler(
    filename = function() {
      if(input$checkGroup2=="all"){
        paste("lasso_EM",input$checkGroup2 ,"results.rds",sep="")}else{
          paste("lasso_EM",input$checkGroup2 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup2=="all"){
        saveRDS(result(), file)}else if(input$checkGroup2=="item"){
          #m<-cbind(result()$ra,result()$rb)
          #colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
          m<-cbind(result()$Beta)
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup2=="cov"){
          #m<-result()$rsigma
          #rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
          m<-result()$Covs
          write.csv(m,file,row.names = F)
        }
    }
  )
  
  ##########################################################################
  ###########################################################################
  ###########################################################################
  #Lasso-EMM output
  resp1<-reactive({req(input$file1)
    df <- data.matrix(read.csv(input$file1$datapath,
                               header = input$header1,
                               sep = input$sep1))
    return(df)
  })
  output$contents1 <- renderTable({
    if(input$disp1 == "head") {
      return(resp1()[1:6,])
    }
    else {
      return(resp1())
    }
    
  })
  ##############################################################
  N.vec1<-reactive({req(input$file2)
    df2 <- data.matrix(read.csv(input$file2$datapath,
                                header = input$header2,
                                sep = input$sep2))
    return(df2)
  })
  
  output$contents2 <- renderTable({
    
    return(N.vec1())
    
  })
  ##############################################################
  m1<-reactive({req(input$m)})
  
  r1<-reactive({req(input$r)})
  y1<-reactive({req(input$y)})
  eta1<-reactive({req(input$eta)})
  ##############################################################
  gra001<-reactive({req(input$file3)
    df3 <- data.matrix(read.csv(input$file3$datapath,
                                header = input$header3,
                                sep = input$sep3))
    return(df3)
  })
  output$contents3 <- renderTable({
    if(input$disp3 == "head") {
      return(gra001()[1:6,])
    }
    else {
      return(gra001())
    }
    
  })
  ##############################################################
  grd001<-reactive({req(input$file4)
    df4 <- data.matrix(read.csv(input$file4$datapath,
                                header = input$header4,
                                sep = input$sep4))
    return(df4)
  })
  output$contents4 <- renderTable({
    if(input$disp4 == "head") {
      return(grd001()[1:6,])
    }
    else {
      return(grd001())
    }
    
  })
  
  ##############################################################
  grgamma001<-reactive({req(input$file5)
    df5 <- data.matrix(read.csv(input$file5$datapath,
                                header = input$header5,
                                sep = input$sep5))
    return(df5)
  })
  output$contents5 <- renderTable({
    if(input$disp5 == "head") {
      return(grgamma001()[1:6,])
    }
    else {
      return(grgamma001())
    }
    
  })
  
  ##############################################################
  grbeta001<-reactive({req(input$file6)
    df6 <- data.matrix(read.csv(input$file6$datapath,
                                header = input$header6,
                                sep = input$sep6))
    return(df6)
  })
  output$contents6 <- renderTable({
    if(input$disp6 == "head") {
      return(grbeta001()[1:6,])
    }
    else {
      return(grbeta001())
    }
    
  })
  ##############################################################
  Mu.list1<-reactive({req(input$file7)
    df7 <- data.matrix(read.csv(input$file7$datapath,
                                header = input$header7,
                                sep = input$sep7))
    return(df7)
  })
  output$contents7 <- renderTable({
    return( Mu.list1())
  })
  ##############################################################
  Sig.list1<-reactive({req(input$file8)
    df8 <- data.matrix(read.csv(input$file8$datapath,
                                header = input$header8,
                                sep = input$sep8))
    return(df8)
  })
  output$contents8 <- renderTable({
    return(Sig.list1())
  })
  ##############################################################
  result<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    result<-Reg_EMM_DIF(resp1(),m1(),r1(),y1(),N.vec1(),eta1(),1e-3,1e-7,gra001(),grd001(),grbeta001(),grgamma001(),Mu.list1(),Sig.list1()) 
    return(result)
  })
  
  output$par2<-renderTable({
    input$go3
    isolate({
      # m<-cbind(result()$ra,result()$rb)
      m<-cbind(result()$Beta)
      #colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
      return(m)
    })
    
  })
  
  output$cov2<-renderTable({
    input$go3
    isolate({
      m<-result()$Covs
      #rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
      return(m)
    })
    
  })
  
  #output$gic2<-renderText({
  # input$go3
  #  isolate({
  #    paste("GIC =", round(result()$GIC,2))
  #  })
  #})
  
  #output$aic2<-renderText({
  #  input$go3
  #  isolate({
  #    paste("AIC =", round(result()$AIC,2))
  #  })
  #})
  
  output$bic2<-renderText({
    input$go3
    isolate({
      paste("BIC =", round(result()$bic,2))
    })
  })
  #Downloadable csv of selected dataset ----
  output$downloadDataall2 <- downloadHandler(
    filename = function() {
      if(input$checkGroup2=="all"){
        paste("lasso_EM",input$checkGroup2 ,"results.rds",sep="")}else{
          paste("lasso_EM",input$checkGroup2 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup2=="all"){
        saveRDS(result(), file)}else if(input$checkGroup2=="item"){
          #m<-cbind(result()$ra,result()$rb)
          #colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
          m<-cbind(result()$Beta)
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup2=="cov"){
          #m<-result()$rsigma
          #rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
          m<-result()$Covs
          write.csv(m,file,row.names = F)
        }
    }
  )
  ###########################################################################
  ###########################################################################
  ###########################################################################
  # Adaptive Lasso-EM output
  resp1<-reactive({req(input$file1)
    df <- data.matrix(read.csv(input$file1$datapath,
                               header = input$header1,
                               sep = input$sep1))
    return(df)
  })
  output$contents1 <- renderTable({
    if(input$disp1 == "head") {
      return(resp1()[1:6,])
    }
    else {
      return(resp1())
    }
    
  })
  ##############################################################
  N.vec1<-reactive({req(input$file2)
    df2 <- data.matrix(read.csv(input$file2$datapath,
                                header = input$header2,
                                sep = input$sep2))
    return(df2)
  })
  
  output$contents2 <- renderTable({
    
    return(N.vec1())
    
  })
  ##############################################################
  m1<-reactive({req(input$m)})
  
  r1<-reactive({req(input$r)})
  y1<-reactive({req(input$y)})
  eta1<-reactive({req(input$eta)})
  ##############################################################
  gra001<-reactive({req(input$file3)
    df3 <- data.matrix(read.csv(input$file3$datapath,
                                header = input$header3,
                                sep = input$sep3))
    return(df3)
  })
  output$contents3 <- renderTable({
    if(input$disp3 == "head") {
      return(gra001()[1:6,])
    }
    else {
      return(gra001())
    }
    
  })
  ##############################################################
  grd001<-reactive({req(input$file4)
    df4 <- data.matrix(read.csv(input$file4$datapath,
                                header = input$header4,
                                sep = input$sep4))
    return(df4)
  })
  output$contents4 <- renderTable({
    if(input$disp4 == "head") {
      return(grd001()[1:6,])
    }
    else {
      return(grd001())
    }
    
  })
  
  ##############################################################
  grgamma001<-reactive({req(input$file5)
    df5 <- data.matrix(read.csv(input$file5$datapath,
                                header = input$header5,
                                sep = input$sep5))
    return(df5)
  })
  output$contents5 <- renderTable({
    if(input$disp5 == "head") {
      return(grgamma001()[1:6,])
    }
    else {
      return(grgamma001())
    }
    
  })
  
  ##############################################################
  grbeta001<-reactive({req(input$file6)
    df6 <- data.matrix(read.csv(input$file6$datapath,
                                header = input$header6,
                                sep = input$sep6))
    return(df6)
  })
  output$contents6 <- renderTable({
    if(input$disp6 == "head") {
      return(grbeta001()[1:6,])
    }
    else {
      return(grbeta001())
    }
    
  })
  ##############################################################
  Mu.list1<-reactive({req(input$file7)
    df7 <- data.matrix(read.csv(input$file7$datapath,
                                header = input$header7,
                                sep = input$sep7))
    return(df7)
  })
  output$contents7 <- renderTable({
    return( Mu.list1())
  })
  ##############################################################
  Sig.list1<-reactive({req(input$file8)
    df8 <- data.matrix(read.csv(input$file8$datapath,
                                header = input$header8,
                                sep = input$sep8))
    return(df8)
  })
  output$contents8 <- renderTable({
    return(Sig.list1())
  })
  ##############################################################
  result<-reactive({
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Iteration times", value = 0)
    # Close the progress when this reactive exits (even if there's an error)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL) {
      if (is.null(value)) {
        value <- progress$getValue()
        value <- value + (progress$getMax() - value) / 10
      }
      progress$set(value = value, detail = detail)
    }
    result<-Reg_Adaptive_DIF(resp1(),m1(),r1(),y1(),N.vec1(),eta1(),1e-3,1e-7,gra001(),grd001(),grbeta001(),grgamma001(),Mu.list1(),Sig.list1()) 
    return(result)
  })
  
  output$par2<-renderTable({
    input$go3
    isolate({
      # m<-cbind(result()$ra,result()$rb)
      m<-cbind(result()$Beta)
      #colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
      return(m)
    })
    
  })
  
  output$cov2<-renderTable({
    input$go3
    isolate({
      m<-result()$Covs
      #rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
      return(m)
    })
    
  })
  
  #output$gic2<-renderText({
  # input$go3
  #  isolate({
  #    paste("GIC =", round(result()$GIC,2))
  #  })
  #})
  
  #output$aic2<-renderText({
  #  input$go3
  #  isolate({
  #    paste("AIC =", round(result()$AIC,2))
  #  })
  #})
  
  output$bic2<-renderText({
    input$go3
    isolate({
      paste("BIC =", round(result()$bic,2))
    })
  })
  #Downloadable csv of selected dataset ----
  output$downloadDataall2 <- downloadHandler(
    filename = function() {
      if(input$checkGroup2=="all"){
        paste("lasso_EM",input$checkGroup2 ,"results.rds",sep="")}else{
          paste("lasso_EM",input$checkGroup2 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup2=="all"){
        saveRDS(result(), file)}else if(input$checkGroup2=="item"){
          #m<-cbind(result()$ra,result()$rb)
          #colnames(m)<-c(paste("a",1:domain1(),sep=""),"b")
          m<-cbind(result()$Beta)
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup2=="cov"){
          #m<-result()$rsigma
          #rownames(m)<-colnames(m)<-paste("dimension",1:domain1(),sep="")
          m<-result()$Covs
          write.csv(m,file,row.names = F)
        }
    }
  )
}
# Run the app ----
shinyApp(ui, server)

