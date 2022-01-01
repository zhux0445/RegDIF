#MBR three DIF functions
library(shiny)
library(psych)
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
library(RcppEigen)
library(abind)

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
    for (int mm=0; mm<(y-1); mm++){
      if (gam(mm,kk)!=0){
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
      if (gam(mm,kk)==0){
        Gamscoall(mm,kk)=0;
      }
    }
  }
  arma::vec Gamsco=Gamscoall.elem(find(gam!=0));
  int len3=0;
  for (int mm=0; mm<(y-1); mm++){
    if (bet(mm)!=0){
      len3++;
    }
  }
  arma::rowvec Betscoall=bet;
  for (int yy=1; yy<y; yy++){
    Betscoall.subvec(yy-1,yy-1) = sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy),0);
  }
  for (int mm=0; mm<(y-1); mm++){
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


Mstep<-'
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
#include <RcppArmadillo.h>
#include <RcppEigen.h>
using namespace arma;
using namespace Rcpp;

//[[Rcpp::export]]
double soft2(double s, double tau){
  arma::vec v1=zeros<arma::vec>(2);
  v1(1)=abs(s)-tau;
  double val=sign(s)* v1.max();
  return(val);
}

//[[Rcpp::export]]
List Mstep(int j, arma::rowvec ng, arma::mat rgk, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma, double maxtol,arma::mat X, arma::mat yallgroup, int y, int G, int r, int m, int eta)
{
  arma::rowvec d=grd.row(j-1);
  arma::rowvec a=gra.row(j-1);
  arma::mat gam=grgamma.slice(j-1);
  arma::rowvec bet=grbeta.row(j-1);
  arma::cube Pstar=zeros<cube>(G,(m-1),y);
  arma::cube Qstar=zeros<cube>(G,(m-1),y);
  arma::cube P=zeros<cube>(G,m,y); 
  
  int miter=0;
  for(miter=0; miter<500; miter++){
    miter++;
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
      for (int mm=0; mm<(y-1); mm++){
        if (gam(mm,kk)!=0){
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
        if (gam(mm,kk)==0){
          Gamscoall(mm,kk)=0;
        }
      }
    }
    arma::vec Gamsco=Gamscoall.elem(find(gam!=0));
    int len3=0;
    for (int mm=0; mm<(y-1); mm++){
      if (bet(mm)!=0){
        len3++;
      }
    }
    arma::rowvec Betscoall=bet;
    for (int yy=1; yy<y; yy++){
      Betscoall.subvec(yy-1,yy-1) = sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy),0);
    }
    for (int mm=0; mm<(y-1); mm++){
      if (bet(mm)==0){
        Betscoall(mm)=0;
      }
    }
    arma::vec Betsco=Betscoall.elem(find(Betscoall!=0));
    arma::rowvec minusgrad = -join_horiz(Dsco,Asco.t(),Gamsco.t(),Betsco.t());
    
    arma::mat FI =zeros<arma::mat>( minusgrad.n_elem, minusgrad.n_elem);
    
    for (int yy=0; yy<y; yy++){
      ( FI.row(0)).subvec(0,0) = ( FI.row(0)).subvec(0,0)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%square(Pstar.slice(yy)%Qstar.slice(yy))%(1/(P.slice(yy)).col(0)+1/(P.slice(yy)).col(1))));
    }
    for (int kk1=1; kk1<(Asco.n_elem+1); kk1++){
      uvec ind1=find(a!=0);
      int ind2=kk1-1;
      int ind3=ind1(ind2);
      for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
        int ind4=kk2-1;
        int ind5=ind1(ind4);
        for (int yy=0; yy<y; yy++){
          ( FI.row(kk1)).subvec(kk2,kk2) = ( FI.row(kk1)).subvec(kk2,kk2)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X.col(ind3))%(X.col(ind5))%sum(square(PQdif.slice(yy))/P.slice(yy),1)));
        }
      }
      for (int yy=0; yy<y; yy++){
        ( FI.row(kk1)).subvec(0,0) = ( FI.row(kk1)).subvec(0,0)+ sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X.col(ind3))%Pstar.slice(yy)%Qstar.slice(yy)%((PQdif.slice(yy)).col(0)/(P.slice(yy)).col(0)-(PQdif.slice(yy)).col(1)/(P.slice(yy)).col(1)));
        FI(0,kk1)=FI(kk1,0);
      }
    }
    if (len2>0){
      for (int kk=(1+Asco.n_elem); kk<(1+Asco.n_elem+len2); kk++){
        uvec ind1=find(gam!=0);
        int ind2=kk-(1+Asco.n_elem);
        int ind3=ind1(ind2);
        int grpnumber=(ind3+1)-floor((ind3+1)/(y-1))*(y-1);
        if (grpnumber==0){
          int grp=y;
          int dimnumber=(ind3+1)/(y-1);
          for (int kk1=(1+Asco.n_elem); kk1<(1+Asco.n_elem+len2); kk1++){
            uvec ind12=find(gam!=0);
            int ind22=kk1-(1+Asco.n_elem);
            int ind32=ind12(ind22);
            int grpnumber2=(ind32+1)-floor((ind32+1)/(y-1))*(y-1);
            if (grpnumber2==0){
              int grp2=y;
              int dimnumber2=(ind32+1)/(y-1);
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            } else {
              int grp2=grpnumber2+1;
              int dimnumber2=floor((ind32+1)/(y-1))+1;
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            }
          }
          for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
            uvec ind13=find(a!=0);
            int ind23=kk2-1;
            int ind33=ind13(ind23);
            ( FI.row(kk)).subvec(kk2,kk2) = (-sum(ng.subvec(((grp)*G),((grp+1)*G-1)).t()%(X.col(ind33))%(X.col(dimnumber-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
            ( FI.row(kk2)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk2,kk2);
          }
          ( FI.row(kk)).subvec(0,0)=sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grp-1)%Qstar.slice(grp-1)%((PQdif.slice(grp-1)).col(0)/(P.slice(grp-1)).col(0)-(PQdif.slice(grp-1)).col(1)/(P.slice(grp-1)).col(1)));
          FI(0,kk)=FI(kk,0);
        } else {
          int grp=grpnumber+1;
          int dimnumber=floor((ind3+1)/(y-1))+1;
          for (int kk1=(1+Asco.n_elem); kk1<(1+Asco.n_elem+len2); kk1++){
            uvec ind12=find(gam!=0);
            int ind22=kk1-(1+Asco.n_elem);
            int ind32=ind12(ind22);
            int grpnumber2=(ind32+1)-floor((ind32+1)/(y-1))*(y-1);
            if (grpnumber2==0){
              int grp2=y;
              int dimnumber2=(ind32+1)/(y-1);
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            } else {
              int grp2=grpnumber2+1;
              int dimnumber2=floor((ind32+1)/(y-1))+1;
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            }
          }
          for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
            uvec ind13=find(a!=0);
            int ind23=kk2-1;
            int ind33=ind13(ind23);
            ( FI.row(kk)).subvec(kk2,kk2) = (-sum(ng.subvec(((grp)*G),((grp+1)*G-1)).t()%(X.col(ind33))%(X.col(dimnumber-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
            ( FI.row(kk2)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk2,kk2);
          }
          ( FI.row(kk)).subvec(0,0)=sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grp-1)%Qstar.slice(grp-1)%((PQdif.slice(grp-1)).col(0)/(P.slice(grp-1)).col(0)-(PQdif.slice(grp-1)).col(1)/(P.slice(grp-1)).col(1)));
          FI(0,kk)=FI(kk,0);
        }
      }
    }
    if (len3>0){
      for (int kk=(1+len+len2); kk<(1+len+len2+len3); kk++){
        uvec ind1=find(bet!=0);
        int ind2=kk-(1+len+len2);
        int ind3=ind1(ind2);
        int grpnumber=ind3;
        ( FI.row(kk)).subvec(kk,kk) = (-sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%square((Pstar.slice(grpnumber+1)).col(0)%(Qstar.slice(grpnumber+1)).col(0))%(1/(P.slice(grpnumber+1)).col(0)+1/(P.slice(grpnumber+1)).col(1))));
        FI(kk,0)=FI(kk,kk);
        FI(0,kk)=FI(kk,kk);
        for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
          uvec ind13=find(a!=0);
          int ind23=kk2-1;
          int ind33=ind13(ind23);
          ( FI.row(kk)).subvec(kk2,kk2) = sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%(X.col(ind33))%Pstar.slice(grpnumber+1)%Qstar.slice(grpnumber+1)%((PQdif.slice(grpnumber+1)).col(0)/(P.slice(grpnumber+1)).col(0)-(PQdif.slice(grpnumber+1)).col(1)/(P.slice(grpnumber+1)).col(1)));
          ( FI.row(kk2)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk2,kk2);
        }
        if (len2>0){
          for (int kk3=(1+Asco.n_elem); kk3<(1+Asco.n_elem+len2); kk3++){
            uvec ind14=find(gam!=0);
            int ind24=kk3-(1+Asco.n_elem);
            int ind34=ind14(ind24);
            int grpnumber2=(ind34+1)-floor((ind34+1)/(y-1))*(y-1);
            if (grpnumber2==0){
              int grp=y;
              int dimnumber=(ind34+1)/(y-1);
              if (grp==(grpnumber+2)){
                ( FI.row(kk)).subvec(kk3,kk3) = sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grpnumber+1)%Qstar.slice(grpnumber+1)%((PQdif.slice(grpnumber+1)).col(0)/(P.slice(grpnumber+1)).col(0)-(PQdif.slice(grpnumber+1)).col(1)/(P.slice(grpnumber+1)).col(1)));
                ( FI.row(kk3)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk3,kk3);
              }
            } else {
              int grp=grpnumber2+1;
              int dimnumber=floor((ind34+1)/(y-1))+1;
              if (grp==(grpnumber+2)){
                ( FI.row(kk)).subvec(kk3,kk3) = sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grpnumber+1)%Qstar.slice(grpnumber+1)%((PQdif.slice(grpnumber+1)).col(0)/(P.slice(grpnumber+1)).col(0)-(PQdif.slice(grpnumber+1)).col(1)/(P.slice(grpnumber+1)).col(1)));
                ( FI.row(kk3)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk3,kk3);
              }
            }
          }
        }
      }
    }
    arma::rowvec add=minusgrad*inv(FI);
    d=d+add(0);
    a(find(a!=0))=a(find(a!=0))+add.subvec(1,len);
    if (len2>0){
      arma::mat gam0=gam;
      gam0(find(gam!=0))=gam(find(gam!=0))+add.subvec(len+1,len+len2).t();
      for (int mm=(len+1); mm<(len+len2+1); mm++){
        int ind1=mm-(len+1);
        uvec ind2=find(gam!=0);
        int ind3=ind2(ind1);
        double sgam=gam0(ind3);
        double taugam=-eta/FI(mm,mm);
        add(mm)=soft2(sgam,taugam)-gam(ind3);
      }
      gam(find(gam!=0))=gam(find(gam!=0))+(add.subvec(len+1,len+len2)).t();
    }
    
    if (len3>0){
      arma::rowvec bet0=bet;
      bet0(find(bet!=0))=bet(find(bet!=0))+add.subvec(len+len2+1,len+len2+len3).t();
      for (int mm=(len+len2+1); mm<(len+len2+len3+1); mm++){
        int ind1=mm-(len+len2+1);
        uvec ind2=find(bet!=0);
        int ind3=ind2(ind1);
        double sbet=bet0(ind3);
        double taubet=-eta/FI(mm,mm);
        add(mm)=soft2(sbet,taubet)-bet(ind3);
      }
      bet(find(bet!=0))=bet(find(bet!=0))+add.subvec(len+len2+1,len+len2+len3).t();
    }
    
    if(sum(abs(add))<maxtol){
      break;
    }
  }
  return List::create(Named("d") = d,Named("a") = a,
                      Named("gam") = gam,Named("bet") = bet);
}


//[[Rcpp::export]]
List Mstepadapt(int j, arma::rowvec ng, arma::mat rgk, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::mat grbeta00, arma::cube grgamma, arma::cube grgamma00, double maxtol,arma::mat X, arma::mat yallgroup, int y, int G, int r, int m, int eta, int lam)
{
  arma::rowvec d=grd.row(j-1);
  arma::rowvec a=gra.row(j-1);
  arma::mat gam=grgamma.slice(j-1);
  arma::mat gammle=grgamma00.slice(j-1);
  arma::rowvec bet=grbeta.row(j-1);
  arma::rowvec betmle=grbeta00.row(j-1);
  arma::cube Pstar=zeros<cube>(G,(m-1),y);
  arma::cube Qstar=zeros<cube>(G,(m-1),y);
  arma::cube P=zeros<cube>(G,m,y); 
  
  int miter=0;
  for(miter=0; miter<500; miter++){
    miter++;
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
      for (int mm=0; mm<(y-1); mm++){
        if (gam(mm,kk)!=0){
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
        if (gam(mm,kk)==0){
          Gamscoall(mm,kk)=0;
        }
      }
    }
    arma::vec Gamsco=Gamscoall.elem(find(gam!=0));
    int len3=0;
    for (int mm=0; mm<(y-1); mm++){
      if (bet(mm)!=0){
        len3++;
      }
    }
    arma::rowvec Betscoall=bet;
    for (int yy=1; yy<y; yy++){
      Betscoall.subvec(yy-1,yy-1) = sum(diff(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy),1,1)%Pstar.slice(yy)%Qstar.slice(yy),0);
    }
    for (int mm=0; mm<(y-1); mm++){
      if (bet(mm)==0){
        Betscoall(mm)=0;
      }
    }
    arma::vec Betsco=Betscoall.elem(find(Betscoall!=0));
    arma::rowvec minusgrad = -join_horiz(Dsco,Asco.t(),Gamsco.t(),Betsco.t());
    
    arma::mat FI =zeros<arma::mat>( minusgrad.n_elem, minusgrad.n_elem);
    
    for (int yy=0; yy<y; yy++){
      ( FI.row(0)).subvec(0,0) = ( FI.row(0)).subvec(0,0)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%square(Pstar.slice(yy)%Qstar.slice(yy))%(1/(P.slice(yy)).col(0)+1/(P.slice(yy)).col(1))));
    }
    for (int kk1=1; kk1<(Asco.n_elem+1); kk1++){
      uvec ind1=find(a!=0);
      int ind2=kk1-1;
      int ind3=ind1(ind2);
      for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
        int ind4=kk2-1;
        int ind5=ind1(ind4);
        for (int yy=0; yy<y; yy++){
          ( FI.row(kk1)).subvec(kk2,kk2) = ( FI.row(kk1)).subvec(kk2,kk2)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X.col(ind3))%(X.col(ind5))%sum(square(PQdif.slice(yy))/P.slice(yy),1)));
        }
      }
      for (int yy=0; yy<y; yy++){
        ( FI.row(kk1)).subvec(0,0) = ( FI.row(kk1)).subvec(0,0)+ sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X.col(ind3))%Pstar.slice(yy)%Qstar.slice(yy)%((PQdif.slice(yy)).col(0)/(P.slice(yy)).col(0)-(PQdif.slice(yy)).col(1)/(P.slice(yy)).col(1)));
        FI(0,kk1)=FI(kk1,0);
      }
    }
    if (len2>0){
      for (int kk=(1+Asco.n_elem); kk<(1+Asco.n_elem+len2); kk++){
        uvec ind1=find(gam!=0);
        int ind2=kk-(1+Asco.n_elem);
        int ind3=ind1(ind2);
        int grpnumber=(ind3+1)-floor((ind3+1)/(y-1))*(y-1);
        if (grpnumber==0){
          int grp=y;
          int dimnumber=(ind3+1)/(y-1);
          for (int kk1=(1+Asco.n_elem); kk1<(1+Asco.n_elem+len2); kk1++){
            uvec ind12=find(gam!=0);
            int ind22=kk1-(1+Asco.n_elem);
            int ind32=ind12(ind22);
            int grpnumber2=(ind32+1)-floor((ind32+1)/(y-1))*(y-1);
            if (grpnumber2==0){
              int grp2=y;
              int dimnumber2=(ind32+1)/(y-1);
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            } else {
              int grp2=grpnumber2+1;
              int dimnumber2=floor((ind32+1)/(y-1))+1;
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            }
          }
          for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
            uvec ind13=find(a!=0);
            int ind23=kk2-1;
            int ind33=ind13(ind23);
            ( FI.row(kk)).subvec(kk2,kk2) = (-sum(ng.subvec(((grp)*G),((grp+1)*G-1)).t()%(X.col(ind33))%(X.col(dimnumber-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
            ( FI.row(kk2)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk2,kk2);
          }
          ( FI.row(kk)).subvec(0,0)=sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grp-1)%Qstar.slice(grp-1)%((PQdif.slice(grp-1)).col(0)/(P.slice(grp-1)).col(0)-(PQdif.slice(grp-1)).col(1)/(P.slice(grp-1)).col(1)));
          FI(0,kk)=FI(kk,0);
        } else {
          int grp=grpnumber+1;
          int dimnumber=floor((ind3+1)/(y-1))+1;
          for (int kk1=(1+Asco.n_elem); kk1<(1+Asco.n_elem+len2); kk1++){
            uvec ind12=find(gam!=0);
            int ind22=kk1-(1+Asco.n_elem);
            int ind32=ind12(ind22);
            int grpnumber2=(ind32+1)-floor((ind32+1)/(y-1))*(y-1);
            if (grpnumber2==0){
              int grp2=y;
              int dimnumber2=(ind32+1)/(y-1);
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            } else {
              int grp2=grpnumber2+1;
              int dimnumber2=floor((ind32+1)/(y-1))+1;
              if (grp==grp2){
                ( FI.row(kk)).subvec(kk1,kk1) = (-sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%(X.col(dimnumber2-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
              }
            }
          }
          for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
            uvec ind13=find(a!=0);
            int ind23=kk2-1;
            int ind33=ind13(ind23);
            ( FI.row(kk)).subvec(kk2,kk2) = (-sum(ng.subvec(((grp)*G),((grp+1)*G-1)).t()%(X.col(ind33))%(X.col(dimnumber-1))%sum(square(PQdif.slice(grp-1))/P.slice(grp-1),1)));
            ( FI.row(kk2)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk2,kk2);
          }
          ( FI.row(kk)).subvec(0,0)=sum(ng.subvec((grp*G),((grp+1)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grp-1)%Qstar.slice(grp-1)%((PQdif.slice(grp-1)).col(0)/(P.slice(grp-1)).col(0)-(PQdif.slice(grp-1)).col(1)/(P.slice(grp-1)).col(1)));
          FI(0,kk)=FI(kk,0);
        }
      }
    }
    if (len3>0){
      for (int kk=(1+len+len2); kk<(1+len+len2+len3); kk++){
        uvec ind1=find(bet!=0);
        int ind2=kk-(1+len+len2);
        int ind3=ind1(ind2);
        int grpnumber=ind3;
        ( FI.row(kk)).subvec(kk,kk) = (-sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%square((Pstar.slice(grpnumber+1)).col(0)%(Qstar.slice(grpnumber+1)).col(0))%(1/(P.slice(grpnumber+1)).col(0)+1/(P.slice(grpnumber+1)).col(1))));
        FI(kk,0)=FI(kk,kk);
        FI(0,kk)=FI(kk,kk);
        for (int kk2=1; kk2<(Asco.n_elem+1); kk2++){
          uvec ind13=find(a!=0);
          int ind23=kk2-1;
          int ind33=ind13(ind23);
          ( FI.row(kk)).subvec(kk2,kk2) = sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%(X.col(ind33))%Pstar.slice(grpnumber+1)%Qstar.slice(grpnumber+1)%((PQdif.slice(grpnumber+1)).col(0)/(P.slice(grpnumber+1)).col(0)-(PQdif.slice(grpnumber+1)).col(1)/(P.slice(grpnumber+1)).col(1)));
          ( FI.row(kk2)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk2,kk2);
        }
        if (len2>0){
          for (int kk3=(1+Asco.n_elem); kk3<(1+Asco.n_elem+len2); kk3++){
            uvec ind14=find(gam!=0);
            int ind24=kk3-(1+Asco.n_elem);
            int ind34=ind14(ind24);
            int grpnumber2=(ind34+1)-floor((ind34+1)/(y-1))*(y-1);
            if (grpnumber2==0){
              int grp=y;
              int dimnumber=(ind34+1)/(y-1);
              if (grp==(grpnumber+2)){
                ( FI.row(kk)).subvec(kk3,kk3) = sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grpnumber+1)%Qstar.slice(grpnumber+1)%((PQdif.slice(grpnumber+1)).col(0)/(P.slice(grpnumber+1)).col(0)-(PQdif.slice(grpnumber+1)).col(1)/(P.slice(grpnumber+1)).col(1)));
                ( FI.row(kk3)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk3,kk3);
              }
            } else {
              int grp=grpnumber2+1;
              int dimnumber=floor((ind34+1)/(y-1))+1;
              if (grp==(grpnumber+2)){
                ( FI.row(kk)).subvec(kk3,kk3) = sum(ng.subvec(((grpnumber+2)*G),((grpnumber+3)*G-1)).t()%(X.col(dimnumber-1))%Pstar.slice(grpnumber+1)%Qstar.slice(grpnumber+1)%((PQdif.slice(grpnumber+1)).col(0)/(P.slice(grpnumber+1)).col(0)-(PQdif.slice(grpnumber+1)).col(1)/(P.slice(grpnumber+1)).col(1)));
                ( FI.row(kk3)).subvec(kk,kk) = ( FI.row(kk)).subvec(kk3,kk3);
              }
            }
          }
        }
      }
    }
    arma::rowvec add=minusgrad*inv(FI);
    d=d+add(0);
    a(find(a!=0))=a(find(a!=0))+add.subvec(1,len);
    if (len2>0){
      arma::mat gam0=gam;
      gam0(find(gam!=0))=gam(find(gam!=0))+add.subvec(len+1,len+len2).t();
      for (int mm=(len+1); mm<(len+len2+1); mm++){
        int ind1=mm-(len+1);
        uvec ind2=find(gam!=0);
        int ind3=ind2(ind1);
        double sgam=gam0(ind3);
        double wgam=pow(abs(gammle(ind3)),lam);
        double taugam=-(eta/wgam)/FI(mm,mm);
        add(mm)=soft2(sgam,taugam)-gam(ind3);
      }
      gam(find(gam!=0))=gam(find(gam!=0))+(add.subvec(len+1,len+len2)).t();
    }
    
    if (len3>0){
      arma::rowvec bet0=bet;
      bet0(find(bet!=0))=bet(find(bet!=0))+add.subvec(len+len2+1,len+len2+len3).t();
      for (int mm=(len+len2+1); mm<(len+len2+len3+1); mm++){
        int ind1=mm-(len+len2+1);
        uvec ind2=find(bet!=0);
        int ind3=ind2(ind1);
        double sbet=bet0(ind3);
        double wgam=pow(abs(betmle(ind3)),lam);
        double taubet=-(eta/wgam)/FI(mm,mm);
        add(mm)=soft2(sbet,taubet)-bet(ind3);
      }
      bet(find(bet!=0))=bet(find(bet!=0))+add.subvec(len+len2+1,len+len2+len3).t();
    }
    
    if(sum(abs(add))<maxtol){
      break;
    }
  }
  return List::create(Named("d") = d,Named("a") = a,
                      Named("gam") = gam,Named("bet") = bet);
}
'
sourceCpp(code=Mstep)

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
      grbeta00=as.matrix(starting5new[,(2*domain+2+(y-2)*domain):(2*domain+1+(y-2)*domain+(y-1))])
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
      #estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
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
    for (rr in 1:r){
      for (nn in 1:(y-1)){
        sparsity1[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:(y-1)){
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
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
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
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(2, iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk0=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G,rgky=rgk[(yy*G+1):((yy+1)*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:(y-1)){
      for(k in 1:r){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:(y-1)){
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

Reg_EMM_DIF <- function(resp,Group,indic,eta,Unif=F,eps =1e-3,max.tol=1e-7,r,y,N.vec,gra00,grd00,grbeta00,grgamma00,Mu.list,Sig.list)
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
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
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
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    for (j in 1:J){
      rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
      #estj=M_step(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grbeta=grbeta,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstep(j=j,ng=ng,rgk=rgk,gra=gra,grd=grd,grbeta=grbeta,grgamma=grgamma,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk0=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G,rgky=rgk[(yy*G+1):((yy+1)*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:(y-1)){
      for(k in 1:r){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:(y-1)){
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
      #estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=lam)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstepadapt(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=eta,lam=lam,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
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
    for (rr in 1:r){
      for (nn in 1:(y-1)){
        sparsity1[nn,rr,j]=ifelse(grgamma[nn,rr,j]==0,0,1)
      }
    }
  }
  sparsity2=grbeta
  for (j in 1:J){
    for (rr in 1:(y-1)){
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
    Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
    Sig.est.slice=array(0,c(r,r,y))
    for (yy in 1:y){
      Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
    }
    #LiA=E_step1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N1=N1,N2=N2,N3=N3,N=N)
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
      #estj=M_step_Adaptive(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,max.tol=max.tol,X=X,y.allgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=lam)
      #gra[j,] <- estj[m:(m+r-1)]*Tau  # re-scale a and gamma
      #grd[j,] <- estj[1:(m-1)]
      #grgamma[,,j] <- matrix(estj[(m+r):(m+r+r*(y-1)-1)],y-1,r)*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      #grbeta[j,] <- estj[(m+r+r*(y-1)):(m+r+r*(y-1)+(m-1)*(y-1)-1)]
      estj=Mstepadapt(j=j,ng=ng,rgk=rgk,grd=grd,gra=gra,grgamma=grgamma,grgamma00=grgamma00,grbeta=grbeta,grbeta00=grbeta00,maxtol=max.tol,X=X,yallgroup=y.allgroup,y=y,G=G,m=m,eta=0,lam=lam,r=r)
      gra[j,] <- estj$a*Tau  # re-scale a and gamma
      grd[j,] <- estj$d
      grgamma[,,j] <- estj$gam*matrix(rep(Tau,(y-1)),y-1,r,byrow = T)
      grbeta[j,] <- estj$bet
    }
    df.d <- abs(dold-grd)
    df.a <- abs(aold-gra)
    df.beta <- abs(betaold-grbeta)
    df.gamma <- abs(gammaold-grgamma)
    iter <- iter+1
    #print(c(2,iter,max(df.a), max(df.d), max(df.beta), max(df.gamma)))
  }
  # AIC BIC
  Mu.est.mat=matrix(Mu.est,y,r,byrow = T)
  Sig.est.slice=array(0,c(r,r,y))
  for (yy in 1:y){
    Sig.est.slice[,,yy]=Sig.est[((yy-1)*r+1):((yy-1)*r+r),]
  }
  LiA=Estep1(resp=resp2,Nvec=N.vec,X=X,y=y,G=G,yallgroup=y.allgroup,Mulist=Mu.est.mat,Siglist=Sig.est.slice,gra=gra, grd=grd, grbeta=grbeta, grgamma=grgamma,r=r,J=J,m=m,N=N)
  ng=ngest(LiA=LiA,y=y,Nvec=N.vec,G=G)
  lh=numeric(J)#likelihood function for each item (overall likelihood by sum over j)
  for (j in 1:J){
    rgk=rgkest(j=j,Xijk=Xijk,LiA=LiA,y=y,Nvec=N.vec,G=G,N=N,m=m)
    sumoverk0=sumoverk(G=G,rgky=rgk[(G+1):(2*G),],aj=gra[j,],dj=grd[j,],betjy=0,gamjy=numeric(r),X=X)
    for (yy in 2:y){
      sumoverk0=sumoverk0+sumoverk(G=G,rgky=rgk[(yy*G+1):((yy+1)*G),],aj=gra[j,],dj=grd[j,],betjy=grbeta[j,(yy-1)],gamjy=grgamma[(yy-1),,j],X=X)
    }
    temp=sumoverk0#-eta*norm(as.matrix(x2), type = "1")##sum over g
    lh[j]=temp
  }
  l0norm1=l0norm2=0
  for(i in 1:J){
    for(j in 1:(y-1)){
      for(k in 1:r){
        l0norm1=l0norm1+(grgamma[j,k,i]!=0)
      }
    }
  }
  for(i in 1:J){
    for(j in 1:(y-1)){
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
  grbetamle=init$grbeta00
  grgammamle=init$grgamma00
  Mu.list=init$Mu0
  Sig.list=init$Sigma0
  
  person=nrow(resp)
  item=ncol(resp)
  domain=nrow(indic)
  y=length(unique(Group)) 
  lbd.center=10+0.033*mean(N.vec)
  lbd.vec=seq((lbd.center-10),(lbd.center+10),4)
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
      sim=Reg_Adaptive_DIF(resp=resp,indic=indic,eta=lbd,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
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
      sim=Reg_Adaptive_DIF(resp=resp,indic=indic,eta=lbd,lam=1,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbetamle,grgamma00=grgammamle,Mu.list=Mu.list,Sig.list=Sig.list)
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
    
    if ((lbd==min(lbd.vec))|(lbd==max(lbd.vec))){
      if (lbd==min(lbd.vec)){
        lbd.vec2=seq(max(0,min(lbd.vec)-12),min(lbd.vec)-4,4)
        bics2=gics2=rep(0,length(lbd.vec2))
        ADmat2=array(double(item*(domain+1)*length(lbd.vec2)),dim = c(item,(domain+1),length(lbd.vec2)))
        Gammas2=array(double((y-1)*domain*item*length(lbd.vec2)),dim = c((y-1),domain,item,length(lbd.vec2)))
        Betas2=array(double(item*(y-1)*length(lbd.vec2)),dim = c(item,(y-1),length(lbd.vec2)))
        Mus2=matrix(0,domain*y,length(lbd.vec2))
        Sigs2=array(double(domain*domain*y*length(lbd.vec2)),dim = c(domain*y,domain,length(lbd.vec2)))
        for (k in 1:length(lbd.vec2))
        {
          if (is.function(updateProgress)) {
            text <- paste0("k=:", k)
            updateProgress(detail = text)
          }
          lbd2=lbd.vec2[k]
          #ptm <- proc.time()
          if (Method=="EM"){
            sim=Reg_DIF(resp=resp,indic=indic,eta=lbd2,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
          } 
          if (Method=="EMM"){
            sim=Reg_EMM_DIF(resp=resp,indic=indic,eta=lbd2,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
          } 
          if (Method=="Adapt"){
            sim=Reg_Adaptive_DIF(resp=resp,indic=indic,eta=lbd2,lam=1,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbetamle,grgamma00=grgammamle,Mu.list=Mu.list,Sig.list=Sig.list)
          }
          #print(proc.time() - ptm)
          bics2[k]=sim$bic
          Gammas2[,,,k]=sim$Gamma
          Betas2[,,k]=sim$Beta
          ADmat2[,,k]=sim$est
          Mus2[,k]=sim$means
          Sigs2[,,k]=sim$Covs
        }
      }
      
      if (lbd==max(lbd.vec)){
        lbd.vec2=seq(max(lbd.vec)+4,max(lbd.vec)+12,4)
        bics2=gics2=rep(0,length(lbd.vec2))
        ADmat2=array(double(item*(domain+1)*length(lbd.vec2)),dim = c(item,(domain+1),length(lbd.vec2)))
        Gammas2=array(double((y-1)*domain*item*length(lbd.vec2)),dim = c((y-1),domain,item,length(lbd.vec2)))
        Betas2=array(double(item*(y-1)*length(lbd.vec2)),dim = c(item,(y-1),length(lbd.vec2)))
        Mus2=matrix(0,domain*y,length(lbd.vec2))
        Sigs2=array(double(domain*domain*y*length(lbd.vec2)),dim = c(domain*y,domain,length(lbd.vec2)))
        for (k in 1:length(lbd.vec2))
        {
          if (is.function(updateProgress)) {
            text <- paste0("k=:", k)
            updateProgress(detail = text)
          }
          lbd2=lbd.vec2[k]
          #ptm <- proc.time()
          if (Method=="EM"){
            sim=Reg_DIF(resp=resp,indic=indic,eta=lbd2,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
          } 
          if (Method=="EMM"){
            sim=Reg_EMM_DIF(resp=resp,indic=indic,eta=lbd2,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbeta00,grgamma00=grgamma00,Mu.list=Mu.list,Sig.list=Sig.list)
          } 
          if (Method=="Adapt"){
            sim=Reg_Adaptive_DIF(resp=resp,indic=indic,eta=lbd2,lam=1,Group=Group,Unif=Unif,r=r,y=y,N.vec=N.vec,gra00=gra00,grd00=grd00,grbeta00=grbetamle,grgamma00=grgammamle,Mu.list=Mu.list,Sig.list=Sig.list)
          }
          #print(proc.time() - ptm)
          bics2[k]=sim$bic
          Gammas2[,,,k]=sim$Gamma
          Betas2[,,k]=sim$Beta
          ADmat2[,,k]=sim$est
          Mus2[,k]=sim$means
          Sigs2[,,k]=sim$Covs
        }
      }
      kk=which.min(c(bics,bics2))
      lbd=c(lbd.vec,lbd.vec2)[kk]
      lbd.vec=c(lbd.vec,lbd.vec2)
      Gammas=abind(Gammas,Gammas2,along = 4)
      Betas=abind(Betas,Betas2,along =3)
      ADmat=abind(ADmat,ADmat2,along =3)
      Mus=cbind(Mus,Mus2)
      Sigs=abind(Sigs,Sigs2,along =3)
      bics=c(bics,bics2)
    }
    #Gamma=array(double((y-1)*domain*item),dim = c((y-1),domain,item,1))
    Gamma=Gammas[,,,kk,drop=F]
    #Gamma=Gammas[,,,kk]
    Beta=Betas[,,kk]
    Amat=ADmat[,1:domain,kk]
    Dmat=ADmat[,domain+1,kk]
    Mu=Mus[,kk]
    Sig=Sigs[,,kk]
    ICs=bics
    IC=bics[kk]
  
  return(list(lbd=lbd,lbd.vec=lbd.vec,Gamma=Gamma,Beta=Beta,Amat=Amat,Dmat=Dmat,Mu=Mu,Sig=Sig,ICs=ICs,IC=IC,domain=domain,y=y))
}



LRT_function=function(resp,indic,Group,Unif=T,updateProgress=NULL){
  if (min(resp)==0){
    resp2=as.matrix(resp)
    resp=resp+1
  } else {
    resp2=as.matrix(resp)-1
  }
  resp=resp2

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

    #md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('slopes'))
  rownames(indic)=paste0("a",1:domain)
  if (Unif==T){
    anchors0=c(1:J)
    diff=1
    while(diff>0){
      md.cons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors0]))
      d=DIF(md.cons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'drop')
      anchors=which(d$adj_pvals>0.05)
      diff=length(anchors0)-length(anchors)
      anchors0=anchors
      #anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
    }
    md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[anchors]))
    dif1=DIF(md.noncons0, which.par = c('d'), p.adjust = 'fdr',scheme = 'add',items2test=c(1:J)[-anchors])
    dif1.t=dif1[which(dif1$adj_pvals<0.05),]
    
    #refit
    if (length(rownames(dif1.t))==0){
      md.refit02 <-multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[1:J]))
    } else {
      md.refit02 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var','slopes',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
    }
    #Gamma=Gamma,Beta=Beta,Amat=Amat,Dmat=Dmat,Mu=Mu,Sig=Sig,domain=domain,y=y
    est=cbind(coef(md.refit02,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      est=cbind(est,coef(md.refit02,simplify=T)[[yy]]$items[,domain+1]-coef(md.refit02,simplify=T)[[1]]$items[,domain+1])
    }
    gra.est=as.matrix(est[,1:domain])
    grd.est=matrix(est[,(domain+1)],J,1)
    grgamma.est=array(0,dim=c((y-1),domain,J))
    grbeta.est=as.matrix(est[,(domain+1+1):(domain+1+1+y-1-1)])
    Sigma.est=matrix(0,domain*y,domain)
    Mu.est = numeric(domain*y)
    for (yy in 1:y){
      Sigma.est[((yy-1)*domain+1):(yy*domain),]=coef(md.refit02,simplify=T)[[yy]]$cov
      Mu.est[((yy-1)*domain+1):(yy*domain)]=coef(md.refit02,simplify=T)[[yy]]$means
    }
  } else {
    anchors0=c(1:J)
    diff=1
    while(diff>0){
      md.cons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors0]))
      d=DIF(md.cons0, which.par = c(rownames(indic),'d'), p.adjust = 'fdr',scheme = 'drop')
      anchors=which(d$adj_pvals>0.05)
      diff=length(anchors0)-length(anchors)
      anchors0=anchors
      #anchors=c(1:20)[-which(colnames(resp)%in%rownames(d))]
    }
    md.noncons0 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[anchors]))
    Jt=c(1:J)[-anchors]
    dif1.t=NULL
    for (jj in Jt){
      an=rownames(indic)[which(indic[,jj]==1)]
      dif1=DIF(md.noncons0, which.par = c(an,'d'), p.adjust = 'fdr',scheme = 'add', items2test=jj)
      dif1.t=rbind(dif1.t,dif1[which(dif1$adj_pvals<0.05),])
    }
    
    #refit
    if (length(rownames(dif1.t))==0){
      md.refit02 <-multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[1:J]))
    } else {
      md.refit02 <- multipleGroup(resp, model, group = Group,SE=TRUE,invariance=c('free_means', 'free_var',colnames(resp)[which(colnames(resp)%in%rownames(dif1.t)==0)]))
    }
    
    est=cbind(coef(md.refit02,simplify=T)[[1]]$items[,1:(domain+m-1)])
    for (yy in 2:y){
      est=cbind(est,coef(md.refit02,simplify=T)[[yy]]$items[,1:domain]-coef(md.refit02,simplify=T)[[1]]$items[,1:domain])
    }
    for (yy in 2:y){
      est=cbind(est,coef(md.refit02,simplify=T)[[yy]]$items[,domain+1]-coef(md.refit02,simplify=T)[[1]]$items[,(domain+1)])
    }
    gra.est=as.matrix(est[,1:domain])
    grd.est=matrix(est[,domain+1],J,1)
    grgamma.est=array(double((y-1)*domain*J),dim = c((y-1),domain,J))
    for (yy in 1:(y-1)){
      grgamma.est[yy,,]=t(est[,(domain+1+(yy-1)*domain+1):(domain+1+(yy-1)*domain+domain)])
    }
    grbeta.est=as.matrix(est[,(domain+1+(y-1)*domain+1):(domain+1+(y-1)*domain+1+y-1-1)])
    Sigma.est=matrix(0,domain*y,domain)
    Mu.est = numeric(domain*y)
    for (yy in 1:y){
      Sigma.est[((yy-1)*domain+1):(yy*domain),]=coef(md.refit02,simplify=T)[[yy]]$cov
      Mu.est[((yy-1)*domain+1):(yy*domain)]=coef(md.refit02,simplify=T)[[yy]]$means
    }
  }
    #gra.est=matrix(0,J,domain)
    #grd.est=matrix(0,J,1)
    #grgamma.est=array(0,dim=c((y-1),domain,J))
    #grbeta.est=matrix(0,J,y-1)
    #Sigma.est=matrix(0,domain*y,domain)
    #Mu.est = numeric(domain*y)
  
  return(list(Gamma=grgamma.est,Beta=grbeta.est,Amat=gra.est,Dmat=grd.est,Mu=Mu.est,Sig= Sigma.est,domain=domain,y=y))
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
                                                         "Adaptive lasso EM"="Adapt",
                                                         "Likelihood Ratio Test"="LRT"),
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
                              h3(span(textOutput("warn"),style="color:red")),
                              
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
        value <- value + (progress$getMax() - value) /4
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
    if (input$method == "LRT") {
      result=LRT_function(u(),indic(),Group(),Unif)
    } else {
      result=reg_DIF_alllbd(u(),indic(),Group(),Method,Unif,updateProgress)
    }
    return(result)
  })
  
  output$warn<-renderText({
    input$go1
    isolate({
      if (input$method == "LRT") {
        return(NULL)
      } else {
        if(result0()$lbd == min(result0()$lbd.vec) || result0()$lbd == max(result0()$lbd.vec)){
          return("Warning: The optimal penalty parameter may be out of range")
        }else{
          return(NULL)
        }
      }
      })
  })
  #(lbd=lbd,Gamma=Gamma,Beta=Beta,Amat=Amat,Dmat=Dmat,Sig=Sig,ICs=ICs)
  output$par1<-renderTable({
    input$go1
    isolate({
      if (input$method == "LRT"){
        domain=result0()$domain
        y=result0()$y
        m<-cbind(result0()$Amat,result0()$Dmat)
          for (r in 1:domain){
            m<-cbind(m,t(result0()$Gamma[,r,]))
          }
        m<-cbind(m,result0()$Beta)
        gp=NULL
        #for(yy in 1:(y-1)){
        #  gp=c(gp,rep(yy,domain))
        #}
        #colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))
        #rownames(m)<-c(paste("Item",1:ncol(indic),sep=""))
        for(r in 1:domain){
          gp=c(gp,rep(r,y-1))
        }
        colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:(y-1),sep=""),result0()$domain),gp,sep=""),paste("beta",1:(y-1),sep=""))
        
      } else {
        domain=result0()$domain
        y=result0()$y
        m<-cbind(result0()$Amat,result0()$Dmat)
        if (y==2){
          for (r in 1:domain){
            m<-cbind(m,result0()$Gamma[,r,,1])
          }
        } else {
          for (r in 1:domain){
            m<-cbind(m,t(result0()$Gamma[,r,,1]))
          }
        }
        m<-cbind(m,result0()$Beta)
        gp=NULL
        #for(yy in 1:(y-1)){
        #  gp=c(gp,rep(yy,domain))
        #}
        #colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))
        #rownames(m)<-c(paste("Item",1:ncol(indic),sep=""))
        for(r in 1:domain){
          gp=c(gp,rep(r,y-1))
        }
        colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:(y-1),sep=""),result0()$domain),gp,sep=""),paste("beta",1:(y-1),sep=""))
      }
      return(m)
    })
    
  }, rownames = T)
  
  output$mean1<-renderTable({
    input$go1
    isolate({
      m<-matrix(result0()$Mu)
      rownames(m)<-paste("group",rep(1:(result0()$y),each=result0()$domain)," dimension",rep(1:(result0()$domain),result0()$y),sep="")
      return(m)
    })
    
  }, rownames = T)
   
  output$cov1<-renderTable({
    input$go1
    isolate({
      m<-result0()$Sig
      colnames(m)<-paste("dimension",1:(result0()$domain),sep="")
      rownames(m)<-paste("group",rep(1:(result0()$y),each=result0()$domain)," dimension",rep(1:(result0()$domain),result0()$y),sep="")
      return(m)
    })
    
  }, rownames = T)
  
  output$plot <- renderPlot({
    input$go1
    isolate({
      if (input$method == "LRT"){
        return(NULL)
      } else {
        bic=result0()$ICs
        eta=result0()$lbd.vec
        plot(eta,bic)
      }
    })
  })
  #Downloadable csv of selected dataset ----
  
  output$downloadData <- downloadHandler(
    filename = function() {
      if(input$checkGroup1=="all"){
        paste("RegDIF",input$checkGroup1 ,"results.rds",sep="")}else{
          paste("RegDIF",input$checkGroup1 ,"results.csv",sep="")}
      
    },
    content = function(file) {
      if(input$checkGroup1=="all"){
        saveRDS(result0(), file)
        }else if(input$checkGroup1=="item"){
          domain=result0()$domain
          y=result0()$y
          m<-cbind(result0()$Amat,result0()$Dmat)
          if (y==2){
            for (r in 1:domain){
              m<-cbind(m,result0()$Gamma[,r,,1])
            }
          } else {
            for (r in 1:domain){
              m<-cbind(m,t(result0()$Gamma[,r,,1]))
            }
          }
          m<-cbind(m,result0()$Beta)
          gp=NULL
          #for(yy in 1:(y-1)){
          #  gp=c(gp,rep(yy,domain))
          #}
          #colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:domain,sep=""),(y-1)),gp,sep=""),paste("beta",1:(y-1),sep=""))
          for(r in 1:domain){
            gp=c(gp,rep(r,y-1))
          }
          colnames(m)<-c(paste("a",1:(result0()$domain),sep=""),"b",paste(rep(paste("gamma",1:(y-1),sep=""),result0()$domain),gp,sep=""),paste("beta",1:(y-1),sep=""))
          rownames(m)<-c(paste("Item",1:ncol(indic()),sep=""))
          write.csv(m,file)
        }else if(input$checkGroup1=="cov"){
          m1<-matrix(result0()$Mu)
          m2<-result0()$Sig
          m=cbind(m1,m2)
          colnames(m)<-c(paste("Mean"),paste("Var dimension",1:(result0()$domain),sep=""))
          rownames(m)<-paste("group",rep(1:(result0()$y),each=result0()$domain)," dimension",rep(1:(result0()$domain),result0()$y),sep="")
          write.csv(m,file)
        }
    }
  )
}
# Run the app ----
shinyApp(ui, server)


