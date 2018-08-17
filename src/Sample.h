#ifndef SAMPLE_H
#define SAMPLE_H


#include <RcppArmadillo.h> 
#include <iostream>
#include <vector>


using namespace Rcpp;

int sample_int(const IntegerVector &x, arma::rowvec prob);
int sample_int_fixed_prob(const IntegerVector &x);
arma::uvec sample_with_replacement_fixed_prob(const IntegerVector &x, const int &N);

#endif
