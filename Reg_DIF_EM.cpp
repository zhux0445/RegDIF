#include<RcppArmadillo.h>

using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double sumoverk(int G, arma::mat rgky, arma::rowvec aj, arma::rowvec dj, arma::rowvec gamjy, arma::mat X){
  double sumoverky=0;
  for (int g =0; g < G; g++ ){
    arma::mat z(1,2,fill::zeros);
    double z_1 = 1-1/(1+exp(-(dj+aj*X[g]+gamjy*X[g])))
    z(1,1) = log(z1);
    z(1,2) = log(1/(1+exp(-(dj+aj*X[g]+gamjy*X[g]))));
    sumoverky += (rgky[g]*z);
  }
  return (sumoverky);
}


