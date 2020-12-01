#include<RcppArmadillo.h>
#include <RcppEigen.h>

using namespace arma;

// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]

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
arma:: vec ngest (arma::mat LiA, int y, arma:: vec Nvec, int G){
  arma:: mat Pi = sum(LiA,1);
  arma:: mat Pirep = repelem(Pi, 1, G);
  arma:: vec ng = sum(LiA/Pirep,0);
  arma:: vec ngallgrp = zeros<arma::vec>(G*y);
  for (int yy =0; yy < y; yy++ ){
    ngallgrp.subvec(((yy-1)*G),((yy-1)*G+G-1))=sum(LiA.submat((sum(Nvec.subvec(0,yy-1))-Nvec(yy-1)),0,(sum(Nvec[1:yy])),(G-1))/Pirep,0);
  }
  return(ng.t());
}

ng.allgrp[((yy-1)*G+1):((yy-1)*G+G)]=apply(LiA[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy])),]/Pi[(sum(N.vec[1:yy])-N.vec[yy]+1):(sum(N.vec[1:yy]))],2,sum)

// [[Rcpp::export]]
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}
