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
SEXP eigenMapMatMult(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B){
  Eigen::MatrixXd C = A * B;
  
  return Rcpp::wrap(C);
}
