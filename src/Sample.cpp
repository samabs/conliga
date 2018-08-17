#include "Sample.h"     

using namespace Rcpp;

// basic functionality borrowed from RcppArmadillo
int sample_int(const IntegerVector &x, arma::rowvec prob)
{
    int nOrig = x.size() - 1;
    int i;
    arma::uvec perm = arma::sort_index(prob, 1); //descending sort of index
    prob = arma::sort(prob, 1);  // descending sort of prob
    // cumulative probabilities
    prob = arma::cumsum(prob);
    // compute the sample 
    double rU = unif_rand();
    for (i = 0; i < nOrig; ++i) {
      if (rU <= prob[i]) return x[perm[i]];
    }
    return x[perm[i]];
}

int sample_int_fixed_prob(const IntegerVector &x)
{
  return x.size() * unif_rand();
}

arma::uvec sample_with_replacement_fixed_prob(const IntegerVector &x, const int &N)
{
    arma::uvec samples(N);
    for (int n = 0; n < N; ++n) {
        samples[n] = sample_int_fixed_prob(x);
    }
    return samples;
}
