#include<RcppArmadillo.h>
#include<omp.h>
using namespace arma;

const double log2pi = std::log(2.0*M_PI);
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends("RcppArmadillo")]]

//[[Rcpp::export]]
double loglike(arma::vec p, arma::mat Z){
  double s = 0;
  int i;
  int n = Z.n_rows;
#pragma omp parallel for reduction(+: s) private(i) shared(p, Z, n) schedule(guided)
  for (i = 0; i < n; i++){
    double si = log(log_sum_exp(p%Z.col(i)));
    s += si;
  }
  return(s);
}
