#ifndef INFERENCE_CPP
#define INFERENCE_CPP

#include "StickyHMM.h"
#include "BetaBinom.h"
using namespace Rcpp;



// [[Rcpp::export]]
void inferBBParamsWithPriors(arma::umat counts,
                             int iterations,
                             double means_sigma,
                             double precision_sigma,
                             arma::mat mean_priors,
                             arma::vec precision_prior,
                             std::string run_ref) {
  IndependentBBs ibbs(counts);
  ibbs.inferParams(iterations,
                   means_sigma,
                   precision_sigma,
                   mean_priors,
                   precision_prior,
                   run_ref);
}


// [[Rcpp::export]]
void runStickyHMM(arma::uvec counts,
                  arma::umat chr_mat,
                  int max_states,
                  arma::vec means,
                  double precision,
                  int iterations,
                  arma::rowvec jump_size_prob,
                  double gamma_a,
                  double gamma_scale,
                  std::string sample_ref,
                  int thin,
                  double norm_sigma_s,
                  double norm_sigma_l,
                  double precision_sigma_s,
                  double precision_sigma_l,
                  double hpp_g_a,
                  double hpp_g_b,
                  double hpp_ak_a,
                  double hpp_ak_b,
                  double hpp_r_c,
                  double hpp_r_d,
                  double precision_p_shape,
                  double precision_p_scale,
                  double gg, // gamma if passed as fixed value
                  double rr, // rho if passed as fixed value
                  double apk, // alpha + kappa if passed as fixed value
                  arma::uvec init_hidden_states,
                  arma::vec init_state_rcn) 
{
  
  StickyHMM new_hmm(max_states,
                    counts,
                    chr_mat,
                    means,
                    precision,
                    init_hidden_states,
                    init_state_rcn);
  
  new_hmm.runMCMC(iterations=iterations,
                  jump_size_prob=jump_size_prob,
                  sample_ref=sample_ref,
                  thin=thin,
                  gamma_a=gamma_a, // prior on RCN
                  gamma_scale=gamma_scale, // prior on RCN
                  norm_sigma_s=norm_sigma_s,
                  norm_sigma_l=norm_sigma_l,
                  precision_sigma_s=precision_sigma_s,
                  precision_sigma_l=precision_sigma_l,
                  hpp_g_a=hpp_g_a, // 
                  hpp_g_b=hpp_g_b,
                  hpp_ak_a=hpp_ak_a,
                  hpp_ak_b=hpp_ak_b,
                  hpp_r_c=hpp_r_c,
                  hpp_r_d=hpp_r_d,
                  precision_p_shape=precision_p_shape,
                  precision_p_scale=precision_p_scale,
                  gg=gg,
                  rr=rr,
                  apk=apk);
  
}

#endif