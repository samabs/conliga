#ifndef DATASIMULATION_CPP
#define DATASIMULATION_CPP

#include "StickyHMM.h"
#include "BetaBinom.h"

using namespace Rcpp;

// [[Rcpp::export]]
arma::umat simulateFromIBBs(arma::vec precs, arma::vec m, arma::urowvec total_draws) {
  IndependentBBs ibbs(precs, m);
  ibbs.simulateCounts(total_draws);
  arma::umat counts = ibbs.getCounts();
  Rcout << ibbs.getPrecision(0.1) << std::endl;
  return counts;
}

// [[Rcpp::export]]
arma::mat getBBQuantiles(arma::vec q, arma::uvec counts, arma::vec means, double precision) {
  IndependentBBs ibbs(counts, precision, means);
  return ibbs.getQuantiles(q);
}

// [[Rcpp::export]]
List simStickyHMM(arma::umat chr_mat, int total_draws, double self_trans,
                  int real_num_states, arma::vec means,
                  double precision, double gamma_a=1.5, double gamma_scale=1/1.5)
{
  StickyHMM hmm(chr_mat, means, precision);
  Rcout << "simulating data." << std::endl;
  hmm.simulateData(real_num_states, total_draws, self_trans, gamma_a, gamma_scale);
  Rcout << "getting data." << std::endl;
  arma::mat trans_mat = hmm.getTransMatrix();
  arma::uvec hidden_states = hmm.getHiddenStates();
  arma::vec lociRCN = hmm.getLociCopyNumber();
  arma::uvec loci_counts = hmm.getCounts();
  arma::vec stateRCN = hmm.getStateCopyNumber();
  arma::rowvec init_dist = hmm.getInitDist();
  double mc_loglik = hmm.getMcll();
  double ibbs_loglik = hmm.getIbbsll();
  Rcout << "returning data." << std::endl;
  return List::create( 
    _["counts"]  = loci_counts, 
    _["hidden_states"]  = hidden_states, 
    _["lociRCN"] = lociRCN,
    _["stateRCN"] = stateRCN,
    _["trans_mat"] = trans_mat,
    _["init_dist"] = init_dist,
    _["mc_loglik"] = mc_loglik,
    _["ibbs_loglik"] = ibbs_loglik
  ) ;
}

// [[Rcpp::export]]
List simStickyHMMGivenStateCopyNum(arma::umat chr_mat, int total_draws, double self_trans,
                                   int real_num_states, arma::vec means,
                                   double precision, arma::vec RCNs)
{
  StickyHMM hmm(chr_mat, means, precision);
  Rcout << "simulating data." << std::endl;
  hmm.simulateDataGivenStateCopyNum(real_num_states, self_trans, RCNs, total_draws);
  Rcout << "getting data." << std::endl;
  arma::mat trans_mat = hmm.getTransMatrix();
  arma::uvec hidden_states = hmm.getHiddenStates();
  arma::vec lociCN = hmm.getLociCopyNumber();
  arma::vec lociRCN = hmm.getLociRelativeCopyNumber(lociCN, means);
  arma::uvec loci_counts = hmm.getCounts();
  arma::vec stateCN = hmm.getStateCopyNumber();
  arma::rowvec init_dist = hmm.getInitDist();
  double mc_loglik = hmm.getMcll();
  double ibbs_loglik = hmm.getIbbsll();
  Rcout << "returning data." << std::endl;
  return List::create( 
    _["counts"]  = loci_counts, 
    _["hidden_states"]  = hidden_states, 
    _["lociCN"] = lociCN,
    _["lociRCN"] = lociRCN,
    _["stateCN"] = stateCN,
    _["trans_mat"] = trans_mat,
    _["init_dist"] = init_dist,
    _["mc_loglik"] = mc_loglik,
    _["ibbs_loglik"] = ibbs_loglik
  ) ;
  
}

// [[Rcpp::export]]
List simCountsFromLociCopyNum(arma::vec lcns, arma::umat chr_mat, arma::vec means,
                              double precision, int total_draws, arma::mat trans_mat,
                              arma::rowvec init_dist, arma::uvec hs)
  
{
  StickyHMM hmm(chr_mat, means, precision);
  Rcout << "simulating data." << std::endl;
  hmm.simulateCountsFromLociCopyNum(lcns, total_draws, trans_mat, init_dist, hs);
  Rcout << "getting data." << std::endl;
  arma::uvec loci_counts = hmm.getCounts();
  arma::vec stateCN = hmm.getStateCopyNumber();
  arma::vec loci_RCN = hmm.getLociRelativeCopyNumber(lcns, means);
  double mc_loglik = hmm.getMcll();
  double ibbs_loglik = hmm.getIbbsll();
  Rcout << "returning data." << std::endl;
  return List::create( 
    _["counts"] = loci_counts,
    _["hidden_states"]  = hs,
    _["lociCN"] = lcns,
    _["lociRCN"] = loci_RCN,
    _["stateCN"] = stateCN,
    _["trans_mat"] = trans_mat,
    _["init_dist"] = init_dist,
    _["mc_loglik"] = mc_loglik,
    _["ibbs_loglik"] = ibbs_loglik
  ) ;
}


#endif