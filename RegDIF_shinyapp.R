#MBR three DIF functions
library(shiny)
library(psych)
library(Rcpp)
library(testit)
library(Matrix)
library(gtools)
library(MASS)
library(mirt) 
library(cacIRT)
library(mvtnorm)
library(graphics)
library(dmutate)
library(Rcpp)
library(RcppParallel)
library(RcppArmadillo)

sumoverk<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
// [[Rcpp::export]]
arma::mat sumoverk (int G, arma::mat rgky, arma::rowvec aj, arma::rowvec dj,arma::rowvec betjy, arma::rowvec gamjy, arma::mat X){
  arma::mat sumoverky = zeros<mat>(1,1);
  arma::mat logdiff = zeros<mat>(2,1);
  for (int g =0; g < G; g++ ){
    logdiff = join_vert(log(1-1/(1+exp(-(dj+aj*X.row(g).t()+betjy+gamjy*X.row(g).t())))),log(1/(1+exp(-(dj+aj*X.row(g).t()+betjy+gamjy*X.row(g).t())))));
    sumoverky += (rgky.row(g))*(logdiff);
  }
  return (sumoverky);
}
'
sourceCpp(code=sumoverk)

ngest<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
// [[Rcpp::export]]
arma::rowvec ngest (arma::mat LiA, int y, arma::uvec Nvec, int G){
  arma::mat Pi = sum(LiA,1);
  arma::mat Pirep = repelem(Pi, 1, G);
  arma::rowvec ng = sum(LiA/Pirep,0);
  arma::rowvec ngallgrp = zeros<arma::rowvec>(G*y);
  for (int yy=0; yy<y; yy++){
    ngallgrp.subvec((yy*G),(yy*G+G-1)) = sum(((LiA.submat((sum(Nvec.subvec(0,yy))-Nvec(yy)),0,(sum(Nvec.subvec(0,yy))-1),(G-1)))/(Pirep.submat((sum(Nvec.subvec(0,yy))-Nvec(yy)),0,(sum(Nvec.subvec(0,yy))-1),(G-1)))),0);
  }
  return(join_horiz(ng,ngallgrp));
}
'
sourceCpp(code=ngest)

rgkest<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
// [[Rcpp::export]]
arma::mat rgkest (int j, arma::cube Xijk, arma::mat LiA, int y, arma::vec Nvec, int G, int N, int m){
  arma::mat Pi = sum(LiA,1);
  arma::cube rLiA = zeros<cube>(N,G,m);
  arma::mat Pirep = repelem(Pi, 1, G);
  
  for (int k=0; k<m; k++){
    arma::mat Xjkmat = ((Xijk.slice(k)).col(j-1));
    arma::mat Xijkrep = repelem(Xjkmat, 1, G);
    rLiA.slice(k)= Xijkrep%LiA/Pirep;
  }
  arma::mat rgk = sum( rLiA,0);
  for (int yy=0; yy<y; yy++){
    int Ny=Nvec(yy);
    arma::mat rgky=sum( rLiA.rows(sum(Nvec.subvec(0,yy))-Ny,sum(Nvec.subvec(0,yy))-1),0);
    rgk=join_cols(rgk,rgky);
  }
  return(rgk);
}
'
sourceCpp(code=rgkest)


eigenMapMatMult<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

'
sourceCpp(code=eigenMapMatMult)




Estep1<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
const double log2pi = std::log(2.0*M_PI);

// [[Rcpp::export]]
double dmvnrm2(arma::rowvec x, arma::rowvec mean, arma::mat sigma, bool logd =false){
  int xdim = x.n_elem;
  double out;
  arma::mat upper = inv(trimatu(chol(sigma))).t();
  double uppersum = sum(log(upper.diag()));
  double constants = -(static_cast<double>(xdim)/2.0)*log2pi;
  arma::vec z = upper*(x - mean).t();
  out = constants - 0.5*sum(z%z) + uppersum;
  if (logd == false){
    out = exp(out);
  }
  return(out);
}

// [[Rcpp::export]]
arma::mat Estep1 (arma::mat resp, arma::vec Nvec, arma::mat X, int y, int G, arma::mat yallgroup, arma::mat Mulist, arma::cube Siglist, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int N)
{
  arma::vec Aallgroups = zeros<arma::vec>(G*y);
  for (int yy = 0; yy < y; yy++){
    for (int ii = 0; ii < G; ii++){
      Aallgroups(yy*G+ii)=dmvnrm2(X.row(ii), Mulist.row(yy), Siglist.slice(yy),FALSE);
    }
  }
  
  arma::mat axmat=gra*X.t(); 
  arma::mat ygamallgroups=zeros<mat>(J*y,r);
  for (int yy = 0; yy < y; yy++){
    for (int j = 0; j < J; j++){
      ygamallgroups.row(yy*J+j)=yallgroup.row(yy)*grgamma.slice(j);
    }
  }
  arma::mat ygamatallgrp=zeros<mat>(J*y,X.n_rows);
  for (int yy = 0; yy < y; yy++){
    ygamatallgrp.rows((yy*J),(yy*J+J-1))=ygamallgroups.rows((yy*J),(yy*J+J-1))*X.t(); 
  }
  arma::mat grbetaallgrp=zeros<mat>(J*y,1); 
  for (int yy = 0; yy < y; yy++){
    grbetaallgrp.rows((yy*J),(yy*J+J-1))=(yallgroup.row(yy)*grbeta.t()).t(); 
  }
  
  arma::mat LiA = zeros<mat>(N,G); 
  for (int yy=0; yy<y; yy++)
  {
    int Ny=Nvec(yy);
    arma::cube pstar1=zeros<cube>(J,(m-1),G);
    arma::cube p1=zeros<cube>(J,m,G);
    for (int g = 0; g < G; g++)
    {
      pstar1.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+grbetaallgrp.rows(yy*J,(yy+1)*J-1)+ygamatallgrp(span(yy*J,(yy+1)*J-1),span(g,g)))));
      p1.slice(g) = (-diff(join_horiz(ones<colvec>(J),(pstar1.slice(g)),zeros<colvec>(J)),1,1));
    }
    arma::cube pij = zeros<cube>(J,G,Ny); 
    for (int j = 0; j < J; j++)
    {
      for (int g = 0; g < G; g++)
      {
        for (int n = 0; n < Ny; n++){
          pij(j,g,n)=p1(j,resp(sum(Nvec.subvec(0,yy))-Ny+n,j),g);
        }
      }
    }
    for (int g = 0; g < G; g++)
    {
      for (int n = 0; n < Ny; n++){
        (LiA.rows(sum(Nvec.subvec(0,yy))-Ny,sum(Nvec.subvec(0,yy))-1))(n,g)=prod((pij.slice(n)).col(g))*Aallgroups(yy*G+g);
      }
    }
  }
  return(LiA);
}
'
sourceCpp(code=Estep1)


scocal<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
// [[Rcpp::export]]
arma::mat scocal(int j, arma::rowvec ng, arma::mat rgk, arma::rowvec a, arma::rowvec d, arma::rowvec bet, arma::mat gam, double maxtol,arma::mat X, arma::mat yallgroup, int y, int G, int r, int m, int eta)
{
  arma::cube Pstar=zeros<cube>(G,(m-1),y);
  arma::cube Qstar=zeros<cube>(G,(m-1),y);
  arma::cube P=zeros<cube>(G,m,y); 
  
  
  for(int yy=0; yy<y; yy++){
    for(int g = 0; g < G; g++){
      (Pstar.slice(yy)).row(g)=1/(1+exp(-(d+a*(X.row(g)).t()+ yallgroup.row(yy)*bet.t()+yallgroup.row(yy)*gam*X.row(g).t())));
      (P.slice(yy)).row(g)=-diff(join_horiz(ones<colvec>(1),(Pstar.slice(yy)).row(g),zeros<colvec>(1)),1,1);
    }
  }
  Qstar=ones<cube>(G,(m-1),y)-Pstar;
  arma::rowvec Dsco=zeros<rowvec>(m-1);
  for (int yy=0; yy<y; yy++){
    Dsco += sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy));
  }
  arma::cube PQdif=zeros<cube>(G,m,y); 
  for (int yy=0; yy<y; yy++){
    PQdif.slice(yy)=-diff(join_horiz(zeros<colvec>(G),Pstar.slice(yy)%Qstar.slice(yy),zeros<colvec>(G)),1,1);
  }
  
  int len=0;
  for (int kk=0; kk<r; kk++){
    if (a(kk)!=0){
      len++;
    }
  }
  
  arma::rowvec Ascoall=zeros<rowvec>(a.n_elem);
  arma::rowvec a01=a;
  for (int kk=0; kk<r; kk++){
    if (a(kk)!=0){
      a01(kk)=1;
    }
  }
  for (int yy=0; yy<y; yy++){
    Ascoall += (sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1).t()*X);
  }
  arma::vec Asco=Ascoall(find(a!=0));
  int len2=0;
  for (int kk=0; kk<r; kk++){
    for (int mm=0; mm<2; mm++){
      if (gam(kk,mm)!=0){
        len2++;
      }
    }
  }
  arma::mat Gamscoall=gam;
  for (int yy=1; yy<y; yy++){
    Gamscoall.row(yy-1)=(sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1)).t()*X;
  }
  for (int kk=0; kk<r; kk++){
    for (int mm=0; mm<(y-1); mm++){
      if (gam(kk,mm)==0){
        Gamscoall(kk,mm)=0;
      }
    }
  }
  arma::vec Gamsco=Gamscoall.elem(find(gam!=0));
  int len3=0;
  for (int mm=0; mm<2; mm++){
    if (bet(mm)!=0){
      len3++;
    }
  }
  arma::rowvec Betscoall=bet;
  for (int yy=1; yy<y; yy++){
    Betscoall.subvec(yy-1,yy-1) = sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy),0);
  }
  for (int mm=0; mm<2; mm++){
    if (bet(mm)==0){
      Betscoall(mm)=0;
    }
  }
  arma::vec Betsco=Betscoall.elem(find(Betscoall!=0));
  arma::rowvec minusgrad = -join_horiz(Dsco,Asco.t(),Gamsco.t(),Betsco.t());
  
  return(minusgrad);
}
'
sourceCpp(code=scocal)

soft=function(s, tau) {
  val=sign(s)*max(c(abs(s) - tau,0))
  return(val) }

M_step=function(j,ng,rgk,grd,gra,grgamma,grbeta,max.tol,X,y.allgroup,y,G,m,eta,r){
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
    
    PQdif=array(double(G*m*y),dim=c(G,m,y))
    for (yy in 1:y){
      PQdif[,,yy] <- -t(apply(cbind(0,Pstar[,,yy]*Qstar[,,yy],0),1, diff))
    }
    
    Dsco=numeric(length(d));Asco=numeric(sum(a!=0));Gamsco=numeric(sum(gam!=0));Betsco=numeric(sum(bet!=0))
    minusgrad <- scocal(j=j,ng=ng,rgk=rgk,d=d,a=a,gam=gam,bet=bet,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,r=r,eta=eta)
    
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
    for (kk1 in (length(Dsco)+1):(length(Dsco)+length(Asco))){
      for (kk2 in (length(Dsco)+1):(length(Dsco)+length(Asco))){
        for (yy in 1:y){
          FI[kk1,kk2] = FI[kk1,kk2]+ (-sum(ng[(yy*G+1):(yy*G+G)]*((X[,(which(ifelse(a==0,0,1)!=0))[(kk1-length(Dsco))]]))*X[,(which(ifelse(a==0,0,1)!=0))[(kk2-length(Dsco))]]*apply(PQdif[,,yy]^2/P[,,yy],1,sum)))
        }
      }
      for (bb in 1:length(Dsco)){
        for (yy in 1:y){
          FI[kk1,bb] <-  FI[kk1,bb]+ sum(ng[(yy*G+1):(yy*G+G)]*X[,(which(ifelse(a==0,0,1)!=0))[(kk1-length(Dsco))]]*(Pstar[,bb,yy]*Qstar[,bb,yy])*(PQdif[,bb,yy]/P[,bb,yy]-PQdif[,bb+1,yy]/P[,bb+1,yy]))
        }
        FI[bb,kk1] <-  FI[kk1,bb]
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
        
        for (kk1 in (length(Dsco)+length(Asco)+1):(length(Dsco)+length(Asco)+length(Gamsco))){
          grp.number2=which(gam!=0)[kk1-(length(Dsco)+length(Asco))]%%(y-1)
          if (grp.number2==0){
            grp2=y
            dim.number2=which(gam!=0)[kk1-(length(Dsco)+length(Asco))]/(y-1)
          } else {
            grp2=grp.number2+1
            dim.number2=which(gam!=0)[kk1-(length(Dsco)+length(Asco))]%/%(y-1)+1
          }
          if (grp==grp2){
            FI[kk,kk1] = -sum(ng[(grp*G+1):(grp*G+G)]*X[,dim.number]*X[,dim.number2]*apply(PQdif[,,grp]^2/P[,,grp],1,sum))
          }
        }  
        
        for (kk2 in (length(Dsco)+1):(length(Dsco)+length(Asco))){
          FI[kk,kk2]= FI[kk2,kk]= -sum(ng[(grp*G+1):(grp*G+G)]*X[,(which(ifelse(a==0,0,1)!=0))[(kk2-length(Dsco))]]*X[,dim.number]*apply(PQdif[,,grp]^2/P[,,grp],1,sum))
        }
        
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
        for (kk2 in (length(Dsco)+1):(length(Dsco)+length(Asco))){
          FI[kk,kk2] =  FI[kk2,kk]= sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*X[,(which(ifelse(a==0,0,1)!=0))[(kk2-length(Dsco))]]*Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1]*(PQdif[,1,grp.number+1]/P[,1,grp.number+1]-PQdif[,2,grp.number+1]/P[,2,grp.number+1]))
        }
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
              FI[kk,kk2] =  FI[kk2,kk]= sum(ng[((grp.number+1)*G+1):((grp.number+2)*G)]*(X[,dim.number])*Pstar[,1,grp.number+1]*Qstar[,1,grp.number+1]*(PQdif[,1,grp.number+1]/P[,1,grp.number+1]-PQdif[,2,grp.number+1]/P[,2,grp.number+1])) #2pl only, not for GRM
            }
          }
        }
      }
    }
    
    
    add <- qr.solve(FI,t(minusgrad))
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
M_step_Adaptive=function(j,ng,rgk,grd,gra,grgamma,grgamma00,grbeta,grbeta00,max.tol,X,y.allgroup,y,G,m,eta,lam,r){
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

  DIF_init=function(resp,Group,indic,Unif){
    m=2 ##fixed, 2pl only
    N=nrow(resp)
    J=ncol(resp)
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
    # defalt for no impact (when using mirt to estimate MLE, fix the mean and variance for all groups)
    COV <- matrix(TRUE,domain,domain); diag(COV)=FALSE
    model <- mirt.model(t(indic), COV=COV) ##
    if (Unif==T){
      md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('slopes'))
      starting5new=cbind(coef(md.noncons0,simplify=T)[[1]]$items[,1:(domain+m-1)])
      for (yy in 2:y){
        starting5new=cbind(starting5new,coef(md.noncons0,simplify=T)[[yy]]$items[,domain+1]-coef(md.noncons0,simplify=T)[[1]]$items[,domain+1])
        
      }
      gra00=as.matrix(starting5new[,1:domain])
      rownames(gra00) <- c()
      colnames(gra00) <- c()
      grd00=matrix(starting5new[,domain+1],J,1)
      grgamma00=array(0,dim=c((y-1),domain,J))
      grbeta00=as.matrix(starting5new[,(domain+1+1):(domain+1+1+y-1-1)])
      rownames(grbeta00) <- c()
      colnames(grbeta00) <- c()
      #Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
      Sigma0=matrix(0,domain*y,domain)
      #Mu0 = matrix(0,domain,y)
      Mu0 = numeric(domain*y)
      for (yy in 1:y){
        Sigma0[((yy-1)*domain+1):(yy*domain),]=coef(md.noncons0,simplify=T)[[yy]]$cov
        #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
      }
    } else {
      md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE)
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
      grbeta00=as.matrix(starting5new[,(2*domain+2+(y-2)*domain):(2*domain+1+(y-2)*domain+m)])
      rownames(grbeta00) <- c()
      colnames(grbeta00) <- c()
      #Sigma0=array(double(domain*domain*y),dim = c(domain,domain,y))
      Sigma0=matrix(0,domain*y,domain)
      Mu0 = matrix(0,domain,y)
      for (yy in 1:y){
        Sigma0[((yy-1)*domain+1):(yy*domain),]=coef(md.noncons0,simplify=T)[[yy]]$cov
        #Mu0[,yy]=coef(md.noncons0,simplify=T)[[yy]]$means
      }
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

Reg_DIF <- function(resp,Group,indic,Unif,eta,eps =1e-3,max.tol=1e-7,r,y,N.vec=N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  m=2 #fixed 2pl
  #m,r,y,N.vec,gra00=NULL,grd00=NULL,grbeta00=NULL,grgamma00=NULL,Mu.list=NULL,Sig.list= NULL
 
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
  grgamma=grgamma00
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
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
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
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
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
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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

Reg_EMM_DIF <- function(resp,Group,indic,eta,Unif=F,eps =1e-3,max.tol=1e-7,r,y,N.vec=N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list)
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
  m=2 #fixed 2pl
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
    #LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
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
      estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
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
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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

Reg_Adaptive_DIF <- function(resp,Group,indic,eta,lam=1,Unif=F,eps =1e-3,max.tol=1e-7,r,y,N.vec=N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list)
{
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  N <- nrow(resp)
  J <- ncol(resp)
  m=2 #fixed 2pl
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
    #LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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
      estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=lam,r=r)
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
    #LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
    LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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
      estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=lam,r=r)
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
  #LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
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




reg_DIF_alllbd=function(resp,indic,Group,Method,Unif=F,updateProgress=NULL){
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  init=DIF_init(resp=resp2,Group=Group,indic=indic,Unif=Unif)
  m=2 #fixed 2pl
  r=init$r
  y=init$y
  N.vec=as.vector(table(Group))
  gra00=init$gra00
  grd00=init$grd00
  grbeta00=init$grbeta00
  grgamma00=init$grgamma00
  Mu.list=init$Mu0
  Sig.list=init$Sigma0
  
  person=nrow(resp)
  item=ncol(resp)
  domain=nrow(indic)
  y=length(unique(Group)) 
  lbd.vec=seq(10,20,5)
  bics=gics=rep(0,length(lbd.vec))
  ADmat=array(double(item*(domain+1)*length(lbd.vec)),dim = c(item,(domain+1),length(lbd.vec)))
  Gammas=array(double((y-1)*domain*item*length(lbd.vec)),dim = c((y-1),domain,item,length(lbd.vec)))
  Betas=array(double(item*(y-1)*length(lbd.vec)),dim = c(item,(y-1),length(lbd.vec)))
  Mus=matrix(0,domain*y,length(lbd.vec))
  Sigs=array(double(domain*domain*y*length(lbd.vec)),dim = c(domain*y,domain,length(lbd.vec)))
  for (k in 1:1)
  {
    if (is.function(updateProgress)) {
      text <- paste0("k=:", k)
      updateProgress(detail = text)
    }
    lbd=lbd.vec[k]
    #ptm <- proc.time()
    if (Method=="EM"){
      sim=Reg_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
    } 
    if (Method=="EMM"){
      sim=Reg_EMM_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
    } 
    if (Method=="Adapt"){
      sim=Reg_EMM_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
    }
    #print(proc.time() - ptm)
    bics[k]=sim$bic
    Gammas[,,,k]=sim$Gamma
    Betas[,,k]=sim$Beta
    ADmat[,,k]=sim$est
    Mus[,k]=sim$means
    Sigs[,,k]=sim$Covs
    
    gra00=sim$est[,1:domain]
    grd00=matrix(sim$est[,domain+1])
    grbeta00=sim$Beta
    grgamma00=sim$Gamma
    Mu.list=sim$means
    Sig.list=sim$Covs
  }
  for (k in 2:length(lbd.vec))
  {
    lbd=lbd.vec[k]
    #ptm <- proc.time()
    if (Method=="EM"){
      sim=Reg_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
    } 
    if (Method=="EMM"){
      sim=Reg_EMM_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
    } 
    if (Method=="Adapt"){
      sim=Reg_EMM_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
    }
    #print(proc.time() - ptm)
    bics[k]=sim$bic
    Gammas[,,,k]=sim$Gamma
    Betas[,,k]=sim$Beta
    ADmat[,,k]=sim$est
    Mus[,k]=sim$means
    Sigs[,,k]=sim$Covs
  }
  
  #est=cbind(gra,grd),Gamma=grgamma,Beta=grbeta,iter=iter,bic=BIC, means=Mu.est,Covs=Sig.est)
    kk=which.min(bics)
    lbd=lbd.vec[kk]
    Gamma=Gammas[,,,kk]
    Beta=Betas[,,kk]
    Amat=ADmat[,1:domain,kk]
    Dmat=ADmat[,domain+1,kk]
    Mu=Mus[,kk]
    Sig=Sigs[,,kk]
    ICs=bics
    IC=bics[kk]
  
  return(list(lbd=lbd,Gamma=Gamma,Beta=Beta,Amat=Amat,Dmat=Dmat,Mu=Mu,Sig=Sig,ICs=ICs,IC=IC,domain=domain,y=y))
}


#####shiny app######
options(shiny.maxRequestSize = 30*1024^2)
# Define UI for data upload app ----
ui <- navbarPage("Regularized DIF",
                 
                 # App title ----
                 tabPanel("EM DIF",
                          
                          # Sidebar layout with input and output definitions ----
                          sidebarLayout(
                            
                            # Sidebar panel for inputs ----
                            sidebarPanel(
                              ############################
                              # Input: Data file u ----
                              ############################
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
                              
                              # Input: Select number of rows to display ----
                              radioButtons("disp1", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              ############################
                              # Input: Data file Group ----
                              ############################
                              fileInput("file2", "Choose group indicator CSV File",
                                        multiple = TRUE,
                                        accept = c("text/csv",
                                                   "text/comma-separated-values,text/plain",
                                                   ".csv")),
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
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                              # Input: Select number of rows to display ----
                              radioButtons("disp2", "Display",
                                           choices = c(Head = "head",
                                                       All = "all"),
                                           selected = "head"),
                              
                              # Horizontal line ----
                              tags$hr(),
                              
                             
                              
                              
                              ############################
                              # Input: Data file indic ----
                              ############################
                              fileInput("file4", "Choose loading indicator CSV file",
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
                              
                              
                              
                              #Input: Select Reg_DIF methods ----
                              selectInput("method", 
                                          label = "Choose a regulariztion algorithm",
                                          choices = list("lasso EM"='EM', 
                                                         "lasso EMM"='EMM',
                                                         "Adaptive lasso EM"="Adapt"),
                                          selected = "EM"),
                              
                            
                              #Input: Select information criteria ----
                              selectInput("Type", 
                                          label = "Choose a DIF type",
                                          choices = list("Uniform"='T', 
                                                         "Non-uniform"='F'),
                                          selected = "T"),
                              
                              
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              actionButton("go1", "Run"),
                              
                              
                              # Horizontal line ----
                              tags$hr(),
                              radioButtons("checkGroup1", "Download Results",
                                           choices = list("All Results" = "all",
                                                          "Item Parameters" = "item",
                                                          "Covariance Matrix" = "cov"),
                                           selected = "all"),
                              # Button
                              downloadButton("downloadData", "Download Results")
                              
                            ),
                            
                            # Main panel for displaying outputs ----
                            mainPanel(
                              
                              uiOutput("tab"),
                              uiOutput("tab2"),
                              
                              h2("Data"),
                              # Output: Data file ----
                              tableOutput("contents1"),
                              
                              h2("Group Indicator"),
                              # Output: Data file ----
                              tableOutput("contents2"),
                              
                              h2("Loading indicator"),
                              # Output: Data file ----
                              tableOutput("contents4"),
                              
                              #warning
                              h1(textOutput("warn")),
                              
                              h2("Item Parameter Results"),
                              tableOutput("par1"),
                              
                              h2("Mean Vector"),
                              tableOutput("mean1"),
                              
                              h2("Covariance Matrix"),
                              tableOutput("cov1"),
                              
                              
                              h2("Information Criteria"),
                              plotOutput("plot")
                              
                            )
                            
                          )
                 ))

# Define server logic to read selected file ----
server <- function(input, output,session) {
  url <- a("Regularized DIF (GVEM algorithms)", href="https://www.google.com/")
  output$tab <- renderUI({
    tagList("Other Functions:", url)
  })
  url2 <- a("DIF - likelihood ratio test", href="https://www.google.com/",)
  output$tab2 <- renderUI({
    tagList(url2)
  })
  # data matrix
  u<-reactive({req(input$file1)
    df <- data.matrix(read.csv(input$file1$datapath,
                               header = input$header1,
                               sep = input$sep1))
    return(df)
  })
  output$contents1 <- renderTable({
    if(input$disp1 == "head") {
      return(u()[1:6,])
    }
    else {
      return(u())
    }
    
  })
  # group indicator 
  Group<-reactive({req(input$file2)
    df <- read.csv(input$file2$datapath,
                   header = input$header2,
                   sep = input$sep2)[,1]
    return(df)
  })
  output$contents2 <- renderTable({
    if(input$disp2 == "head") {
      return(Group()[1:6])
    }
    else {
      return(Group())
    }
    
  })
  
  
  # indic
  indic<-reactive({req(input$file4)
    df <- data.matrix(read.csv(input$file4$datapath,
                               header = input$header4,
                               sep = input$sep4))
    return(df)
  })
  output$contents4 <- renderTable({
    return(indic())
  })
  
  result0<-reactive({
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
    if (input$method == "EM") {
      Method="EM"
    } 
    if (input$method == "EMM") {
      Method="EMM"
    } 
    if (input$method == "Adapt") {
      Method="Adapt"
    } 
    if (input$Type == "T") {
      Unif="T"
    } 
    if (input$Type == "F") {
      Unif="F"
    } 
    result=reg_DIF_alllbd(u(),indic(),Group(),Method,Unif,updateProgress)
    return(result)
  })
  
  output$warn<-renderText({
    input$go1
    isolate({
      if(result0 ()$lbd == 10 || result0 ()$lbd == 20){
        return("Warning: The optimal penalty parameter may be out of range, a different gamma value is suggested")
      }else{
        return(NULL)
      }})
  })
  #(lbd=lbd,Gamma=Gamma,Beta=Beta,Amat=Amat,Dmat=Dmat,Sig=Sig,ICs=ICs)
  output$par1<-renderTable({
    input$go1
    isolate({
      domain=result0()$domain
      y=result0()$y
      m<-cbind(result0()$Amat,result0()$Dmat)
      for (r in 1:domain){
        m<-cbind(m,t(result0()$Gamma[r,,]))
      }
      m<-cbind(m,result0()$Beta)
      gp=NULL
      for(yy in 1:(y-1)){
        gp=c(gp,rep(yy,domain))
      }
      colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))
      return(m)
    })
    
  })
  
  output$mean1<-renderTable({
    input$go1
    isolate({
      m<-matrix(result0()$Mu)
      rownames(m)<-paste("group",rep(1:(result0()$y),each=result0()$domain)," dimension",rep(1:(result0()$domain),result0()$y),sep="")
      return(m)
    })
    
  })
  
  output$cov1<-renderTable({
    input$go1
    isolate({
      m<-result0()$Sig
      colnames(m)<-paste("dimension",1:(result0()$domain),sep="")
      rownames(m)<-paste("group",rep(1:(result0()$y),each=result0()$domain)," dimension",rep(1:(result0()$domain),result0()$y),sep="")
      return(m)
    })
    
  })
  
  
  output$plot <- renderPlot({
    input$go1
    isolate({
    bic=result0()$ICs
    eta=seq(10,20,5)
    plot(eta,bic)
    })
  })
  #Downloadable csv of selected dataset ----
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if(input$checkGroup1=="all"){
        paste("GVEMDIF",input$checkGroup1 ,"results.rds",sep="")}else{
          paste("GVEMDIF",input$checkGroup1 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup1=="all"){
        saveRDS(result0(), file)}else if(input$checkGroup1=="item"){
          domain=result0()$domain
          y=result0()$y
          m<-cbind(result0()$Amat,result0()$Dmat)
          for (r in 1:domain){
            m<-cbind(m,t(result0()$Gamma[r,,]))
          }
          m<-cbind(m,result0()$Beta)
          gp=NULL
          for(yy in 1:(y-1)){
            gp=c(gp,rep(yy,domain))
          }
          colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))
          rownames<-c(paste("Item",1:ncol(indic),sep=""))
          write.csv(m,file,row.names = F)
        }else if(input$checkGroup1=="cov"){
          m1<-matrix(result0()$Mu)
          m2<-result0()$Sig
          m=cbind(m1,m2)
          colnames(m)<-c(paste("Mean"),paste("Var dimension",1:(result0()$domain),sep=""))
          rownames(m)<-paste("group",rep(1:(result0()$y),each=result0()$domain)," dimension",rep(1:(result0()$domain),result0()$y),sep="")
          write.csv(m,file,row.names = F)
        }
    }
  )
}
# Run the app ----
shinyApp(ui, server)


#username: uwpmetrics  (or uwpmetrics@gmail.com)
# password: COE_psychometrics21
