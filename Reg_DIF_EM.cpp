#include<RcppArmadillo.h>
#include <RcppEigen.h>

using namespace arma;



// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

const double log2pi = std::log(2.0*M_PI);

void inplace_tri_mat_mult(arma::rowvec &x, arma::mat const &trimat){
  arma::uword const n = trimat.n_cols;
  
  for(unsigned j = n; j-- > 0;){
    double tmp(0.);
    for(unsigned i = 0; i <= j; ++i)
      tmp += trimat.at(i, j) * x[i];
    x[j] = tmp;
  }
}


//[[Rcpp::export]]
double soft2(double s, double tau){
  arma::vec v1=zeros<arma::vec>(2);
  v1(1)=abs(s)-tau;
  double val=sign(s)* v1.max();
  return(val);
}


//[[Rcpp::export]]
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
arma::mat sumoverk (int G, arma::mat rgky, arma::rowvec aj, arma::rowvec dj,arma::rowvec betjy, arma::rowvec gamjy, arma::mat X){
  arma::mat sumoverky = zeros<mat>(1,1);
  arma::mat logdiff = zeros<mat>(2,1);
  for (int g =0; g < G; g++ ){
    logdiff = join_vert(log(1-1/(1+exp(-(dj+aj*X.row(g).t()+betjy+gamjy*X.row(g).t())))),log(1/(1+exp(-(dj+aj*X.row(g).t()+betjy+gamjy*X.row(g).t())))));
    sumoverky += (rgky.row(g))*(logdiff);
  }
  return (sumoverky);
}


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

// [[Rcpp::export]]
arma::mat rgkest0 (int j, arma::cube Xijk, arma::mat LiA, int y, arma::vec Nvec, int G, int N, int m){
  arma::mat Pi = sum(LiA,1);
  arma::cube rLiA = zeros<cube>(N,G,m);
  arma::mat Pirep = repelem(Pi, 1, G);
  
  for (int k=0; k<m; k++){
    arma::mat Xjkmat = ((Xijk.slice(k)).col(j-1));
    arma::mat Xijkrep = repelem(Xjkmat, 1, G);
    rLiA.slice(k)= Xijkrep%LiA/Pirep;
  }
  arma::mat rgk = sum( rLiA,0);
  arma::mat rgky1 = sum( rLiA.rows(0,N/y-1),0);
  arma::mat rgky2 = sum( rLiA.rows(N/y,2*N/y-1),0);
  arma::mat rgky3 = sum( rLiA.rows(2*N/y,3*N/y-1),0);
  return(join_cols(rgk,rgky1,rgky2,rgky3));
}


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


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}


// [[Rcpp::export]]
arma::mat E_step1 (arma::mat resp, arma::vec Nvec, arma::mat X, int y, int G, arma::mat yallgroup, arma::mat Mulist, arma::cube Siglist, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int N)
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


//[[Rcpp::export]]
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

//[[Rcpp::export]]
List M_step(int j, arma::rowvec ng, arma::mat rgk, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma, double maxtol,arma::mat X, arma::mat yallgroup, int y, int G, int r, int m, int eta)
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
    for (int mm=0; mm<(y-1); mm++){
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

