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
double soft(double s, double tau){
  arma::vec v1=zeros<arma::vec>(2);
  v1(1)=abs(s)-tau;
  double val=sign(s)* v1.max();
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
  arma::mat rgky1 = sum( rLiA.rows(0,N/y-1),0);
  arma::mat rgky2 = sum( rLiA.rows(N/y,2*N/y-1),0);
  arma::mat rgky3 = sum( rLiA.rows(2*N/y,3*N/y-1),0);
  return(join_cols(rgk,rgky1,rgky2,rgky3));
}


// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}

// [[Rcpp::export]]
arma::mat E_step1 (arma::mat resp, arma::vec Nvec, arma::mat X, int y, int G, arma::mat yallgroup, arma::mat Mulist, arma::cube Siglist, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int N1, int N2, int N3,int N)
{
  arma::vec Aallgroups = zeros<arma::vec>(G*y);
  for (int yy = 0; yy < y; yy++){
    for (int ii = 0; ii < G; ii++){
      Aallgroups(yy*G+ii)=dmvnrm2(X.row(ii), Mulist.row(yy), Siglist.slice(yy),FALSE);
    }
  }
  
  arma::mat axmat=gra*X.t(); 
  arma::mat ygamallgroups=zeros<mat>(J*y,y-1);
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
      
  arma::cube pstar1=zeros<cube>(J,(m-1),G); 
  arma::cube pstar2=zeros<cube>(J,(m-1),G); 
  arma::cube pstar3=zeros<cube>(J,(m-1),G); 
  arma::cube p1=zeros<cube>(J,m,G); 
  arma::cube p2=zeros<cube>(J,m,G); 
  arma::cube p3=zeros<cube>(J,m,G);
        for (int g = 0; g < G; g++)
        {
          pstar1.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+grbetaallgrp.rows(0,J-1)+ygamatallgrp(span(0,J-1),span(g,g)))));
          p1.slice(g) = (-diff(join_horiz(ones<colvec>(J),(pstar1.slice(g)),zeros<colvec>(J)),1,1));
        }
        for (int g = 0; g < G; g++)
        {
          pstar2.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+grbetaallgrp.rows(J,2*J-1)+ygamatallgrp(span(J,2*J-1),span(g,g)))));
          p2.slice(g) = (-diff(join_horiz(ones<colvec>(J),(pstar2.slice(g)),zeros<colvec>(J)),1,1));
        }
        for (int g = 0; g < G; g++)
        {
          pstar3.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+grbetaallgrp.rows(2*J,3*J-1)+ygamatallgrp(span(2*J,3*J-1),span(g,g)))));
          p3.slice(g) = (-diff(join_horiz(ones<colvec>(J),(pstar3.slice(g)),zeros<colvec>(J)),1,1));
        }
        arma::cube pij = zeros<cube>(J,G,N/3); 
        arma::mat LiA = zeros<mat>(N,G); 
          
          for (int j = 0; j < J; j++)
          {
            for (int g = 0; g < G; g++)
            {
              for (int n = 0; n < N1; n++){
                pij(j,g,n)=p1(j,resp(n,j),g);
              }
            }
          }
          for (int g = 0; g < G; g++)
          {
            for (int n = 0; n < N1; n++){
          (LiA.rows(0,N1-1))(n,g)=prod((pij.slice(n)).col(g))*Aallgroups(g);
            }
          }
          
          for (int j = 0; j < J; j++)
          {
            for (int g = 0; g < G; g++)
            {
              for (int n = 0; n < N2; n++){
                pij(j,g,n)=p2(j,resp(N1+n,j),g);
              }
            }
          }
          for (int g = 0; g < G; g++)
          {
            for (int n = 0; n < N2; n++){
              (LiA.rows(N1,N1+N2-1))(n,g)=prod((pij.slice(n)).col(g))*Aallgroups(G+g);
            }
          }
            
            for (int j = 0; j < J; j++)
            {
              for (int g = 0; g < G; g++)
              {
                for (int n = 0; n < N3; n++){
                  pij(j,g,n)=p3(j,resp(N1+N2+n,j),g);
                }
              }
            }
            for (int g = 0; g < G; g++)
            {
              for (int n = 0; n < N3; n++){
                (LiA.rows(N1+N2,N-1))(n,g)=prod((pij.slice(n)).col(g))*Aallgroups(2*G+g);
              }
            }
  return(LiA);
}


// [[Rcpp::export]]
Rcpp::List M_step(int j, arma::rowvec ng, arma::mat rgk, arma::mat X, int y, int G, arma::mat yallgroup, double maxtol, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int eta)
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
    arma::rowvec Asco=zeros<rowvec>(len);
    arma::rowvec a01=a;
    for (int kk=0; kk<r; kk++){
      if (a(kk)!=0){
        a01(kk)=1;
      }
    }
    for (int yy=0; yy<y; yy++){
      Asco += sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1).t()*(X*a01.t());
    }
    int len2=0;
    for (int kk=0; kk<r; kk++){
      for (int mm=0; mm<2; mm++){
        if (gam(kk,mm)!=0){
          len2++;
        }
      }
    }
    arma::vec Gamsco=zeros<vec>(len2);
    arma::mat Gamscoall=gam;
    for (int yy=1; yy<y; yy++){
      Gamscoall.row(yy-1)=(sum(rgk.rows((yy+1)*G,(yy+2)*G-1)/P.slice(yy)%PQdif.slice(yy),1)).t()*X;
    }
    for (int kk=0; kk<r; kk++){
      for (int mm=0; mm<2; mm++){
        if (gam(kk,mm)==0){
          Gamscoall(kk,mm)=0;
        }
      }
    }
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
    arma::rowvec minusgrad = join_horiz(Dsco,Asco,(Gamscoall*a01.t()).t(),Betscoall);
    arma::rowvec minusgrad2 = minusgrad.elem(find(minusgrad!=0)).t();
    
    arma::mat FI =zeros<arma::mat>( minusgrad.n_elem, minusgrad.n_elem);
    
    for (int yy=0; yy<y; yy++){
      ( FI.row(0)).subvec(0,0) = ( FI.row(0)).subvec(0,0)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%square(Pstar.slice(yy)%Qstar.slice(yy))%(1/(P.slice(yy)).col(0)+1/(P.slice(yy)).col(1))));
    }
    
    for (int yy=0; yy<y; yy++){
      ( FI.row(1)).subvec(1,1) = ( FI.row(1)).subvec(1,1)+ (-sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X*a01.t())%(X*a01.t())%sum(square(PQdif.slice(yy))/P.slice(yy),1)));
      ( FI.row(1)).subvec(0,0) = ( FI.row(1)).subvec(0,0)+ sum(ng.subvec(((yy+1)*G),((yy+2)*G-1)).t()%(X*a01.t())%Pstar.slice(yy)%Qstar.slice(yy)%((PQdif.slice(yy)).col(0)/(P.slice(yy)).col(0)-(PQdif.slice(yy)).col(1)/(P.slice(yy)).col(1)));
      FI(0,1)=FI(1,0);
    }
    
    for (int kk=2; kk<(y+1); kk++){
      ( FI.row(kk)).subvec(kk,kk) = (-sum(ng.subvec((kk*G),((kk+1)*G-1)).t()%(X*a01.t())%(X*a01.t())%sum(square(PQdif.slice(kk-1))/P.slice(kk-1),1)));
      FI(kk,1)=FI(kk,kk);
      FI(1,kk)=FI(kk,kk);
      ( FI.row(kk)).subvec(0,0)=sum(ng.subvec((kk*G),((kk+1)*G-1)).t()%(X*a01.t())%Pstar.slice(kk-1)%Qstar.slice(kk-1)%((PQdif.slice(kk-1)).col(0)/(P.slice(kk-1)).col(0)-(PQdif.slice(kk-1)).col(1)/(P.slice(kk-1)).col(1)));
      FI(0,kk)=FI(kk,0);
    }
    
    for (int kk=(y+1); kk<(2*y); kk++){
      ( FI.row(kk)).subvec(kk,kk) = (-sum(ng.subvec(((kk-2)*G),((kk-1)*G-1)).t()%square(Pstar.slice(kk-3)%Qstar.slice(kk-3))%(1/(P.slice(kk-3)).col(0)+1/(P.slice(kk-3)).col(1))));
      FI(kk,0)=FI(kk,kk);
      FI(0,kk)=FI(kk,kk);
      ( FI.row(kk)).subvec(1,1)=sum(ng.subvec(((kk-2)*G),((kk-1)*G-1)).t()%(X*a01.t())%Pstar.slice(kk-3)%Qstar.slice(kk-3)%((PQdif.slice(kk-3)).col(0)/(P.slice(kk-3)).col(0)-(PQdif.slice(kk-3)).col(1)/(P.slice(kk-3)).col(1)));
      FI(1,kk)=FI(kk,1);
      
    }
   
    arma::rowvec scoall01=ones<arma::rowvec>(2+2*(y-1));
    arma::mat sco01mat = zeros<arma::mat>(2+len2+len3,2+2*(y-1));
    arma::mat scoall01mat(2+2*(y-1),2+2*(y-1),fill::eye);
    arma::rowvec gam01=gam*a01.t();
    gam01.elem(find(gam01!=0)).ones();
    scoall01.subvec(2,y)= gam01;
    arma::rowvec bet01=bet;
    bet01.elem(find(bet01!=0)).ones();
    scoall01.subvec(y+1,2*y-1)= bet01;
    for (int kk=0; kk<(2+len2+len3); kk++){
      int ind1=find(scoall01!=0)(kk);
      sco01mat.row(kk)=scoall01mat.row(ind1);
    }
    arma::mat FI2= zeros<arma::mat>(2+len2+len3,2+len2+len3);
    for (int kk1=0; kk1<(2+len2+len3); kk1++){
      for (int kk2=0; kk2<(2+len2+len3); kk2++){
        FI2(kk1,kk2)= sco01mat.row(kk1)*FI*sco01mat.row(kk2).t();
      }
    }
    
    arma::rowvec add=solve(FI2,minusgrad);
    d=+add.subvec(0,0);
    a.elem(find(a!=0))=+add.subvec(1,1);
    
    arma::rowvec bet=grbeta.row(j-1);
    
    if (len2>0){
      arma::mat gam0=gam;
      gam0.elem(find(gam!=0))=gam.elem(find(gam!=0))+add.subvec(2,2+len2-1);
      for (int mm=2; mm<2+len2; mm++){
        add(mm)=soft(gam0.elem(find(gam!=0))(mm-2),-eta/FI2(mm,mm))-gam.elem(find(gam!=0))(mm-2);
      }
      gam.elem(find(gam!=0))=gam.elem(find(gam!=0))+add.subvec(2,2+len2-1);
    }
    if (len3>0){
      arma::mat bet0=bet;
      bet0.elem(find(bet!=0))=bet.elem(find(bet!=0))+add.subvec(2,2+len3-1);
      for (int mm=2; mm<2+len3; mm++){
        add(mm)=soft(bet0.elem(find(bet!=0))(mm-2),-eta/FI2(mm,mm))-bet.elem(find(bet!=0))(mm-2);
      }
      bet.elem(find(bet!=0))=bet.elem(find(bet!=0))+add.subvec(2,2+len3-1);
    }
    
    if(sum(add)<maxtol){
      break;
    } 
  }
  return Rcpp::List::create(Rcpp::Named("a") = a,
                            Rcpp::Named("d") = d,
                            Rcpp::Named("gam") = gam,
                            Rcpp::Named("bet") = bet);
}
