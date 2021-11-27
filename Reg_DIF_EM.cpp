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


