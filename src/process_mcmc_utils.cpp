#ifndef PROCESS_MCMC_UTILS_CPP
#define PROCESS_MCMC_UTILS_CPP

#include <RcppArmadillo.h>

using namespace Rcpp;


// [[Rcpp::export]]
arma::uvec lociStateDiff(arma::umat &m)
{
  arma::uvec diff(m.n_cols);
  diff.zeros();
  
  for (int iter = 0; iter < m.n_rows; ++iter) {
    for (int locus = 1; locus < m.n_cols; ++locus) {
      if (m(iter,locus) != m(iter,locus-1)) {
        diff[locus] += 1;
      }
    }
  }
  
  return diff;
}


#endif