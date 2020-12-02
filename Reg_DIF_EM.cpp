#include<RcppArmadillo.h>
#include <RcppEigen.h>
#include <omp.h>

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

// [[Rcpp::export]]
arma::vec dmvnrm_arma_mc(arma::mat const &x,  
                         arma::rowvec const &mean,  
                         arma::mat const &sigma, 
                         bool const logd = false,
                         int const cores = 1) {  
  using arma::uword;
  omp_set_num_threads(cores);
  uword const n = x.n_rows, 
    xdim = x.n_cols;
  arma::vec out(n);
  arma::mat const rooti = arma::inv(trimatu(arma::chol(sigma)));
  double const rootisum = arma::sum(log(rooti.diag())), 
    constants = -(double)xdim/2.0 * log2pi, 
    other_terms = rootisum + constants;
  
  arma::rowvec z;
#pragma omp parallel for schedule(static) private(z)
  for (uword i = 0; i < n; i++) {
    z = (x.row(i) - mean);
    inplace_tri_mat_mult(z, rooti);   
    out(i) = other_terms - 0.5 * arma::dot(z, z);     
  }  
  
  if (logd)
    return out;
  return exp(out);
}


// [[Rcpp::export]]
arma::mat sumoverk (int G, arma::mat rgky, arma::rowvec aj, arma::rowvec dj, arma::rowvec gamjy, arma::mat X){
  arma::mat sumoverky = zeros<mat>(1,1);
  arma::mat logdiff = zeros<mat>(2,1);
  for (int g =0; g < G; g++ ){
    logdiff = join_vert(log(1-1/(1+exp(-(dj+aj*X.row(g).t()+gamjy*X.row(g).t())))),log(1/(1+exp(-(dj+aj*X.row(g).t()+gamjy*X.row(g).t())))));
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
arma::mat E_step1 (arma::mat resp, arma::vec Nvec, arma::mat X, int y, int G, arma::mat yallgroup, arma::vec Mulist, arma::mat Siglist, arma::mat gra, arma::mat grd, arma::mat grbeta, arma::cube grgamma,int r, int J, int m, int N1, int N2, int N3,int N){
  arma::vec Aallgroups = zeros<arma::vec>(X.n_rows*y);
  for (int yy = 0; yy < y; yy++){
    Aallgroups.subvec((yy*X.n_rows),(yy*X.n_rows+X.n_rows-1))=dmvnrm_arma_mc(X, Mulist.subvec((yy*r),(yy*r+r-1)), Siglist.rows((yy*r),(yy*r+r-1)),FALSE,4);
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
  arma::cube pstar1=zeros<cube>(J,(m-1),G); 
  arma::cube pstar2=zeros<cube>(J,(m-1),G); 
  arma::cube pstar3=zeros<cube>(J,(m-1),G); 
  arma::cube p1=zeros<cube>(J,m,G); 
  arma::cube p2=zeros<cube>(J,m,G); 
  arma::cube p3=zeros<cube>(J,m,G);
        for (int g = 0; g < G; g++)
        {
          pstar1.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+ygamatallgrp(span(0,J-1),span(g,g)))));
          p1.slice(g) = (-diff(join_vert(ones<rowvec>(J),(pstar1.slice(g)).t(),zeros<rowvec>(J)),0)).t();
        }
        for (int g = 0; g < G; g++)
        {
          pstar2.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+ygamatallgrp(span(J,2*J-1),span(g,g)))));
          p2.slice(g) = (-diff(join_vert(ones<rowvec>(J),(pstar2.slice(g)).t(),zeros<rowvec>(J)),0)).t();
        }
        for (int g = 0; g < G; g++)
        {
          pstar3.slice(g) = 1/(1+exp(-(grd+axmat.col(g)+ygamatallgrp(span(2*J,3*J-1),span(g,g)))));
          p3.slice(g) = (-diff(join_vert(ones<rowvec>(J),(pstar3.slice(g)).t(),zeros<rowvec>(J)),0)).t();
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

