#ifndef STICKY_HMM_H
#define STICKY_HMM_H

#include "BetaBinom.h"
#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <limits>
//#define SSTR( x ) static_cast< std::ostringstream & >( \
//        ( std::ostringstream() << std::dec << x ) ).str()

using namespace Rcpp;

class StickyHMM {
	int max_states;
	int num_chrs;
	int num_loci;
	arma::uvec num_loci_per_chr;
	arma::uvec hidden_states;
	arma::mat trans_prob;
	arma::rowvec init_prob;
	arma::mat trans_count;
	arma::mat M;
	arma::mat M_bar;
	arma::rowvec m_bar_colsums;
	arma::rowvec beta;
	arma::vec w;
	IndependentBBs ibbs;
	arma::vec state_RCNs;
	arma::vec loci_RCNs;
	arma::umat chr_start_end;
	arma::mat loglik_cache;

	// Gamma prior parameters for precision
	double precision_prior_shape;
	double precision_prior_scale;

	// Hyper-parameters
	double alpha;
	double gamma;
	double kappa;
	double rho;
	
	bool gamma_fixed;
	bool rho_fixed;
	bool apk_fixed;

	// Hyper-param priors
	// Gamma(a,b) priors on (alpha + kappa) and gamma
	double hpp_gamma_a;
	double hpp_gamma_b;
	double hpp_aplusk_a;
	double hpp_aplusk_b;
	// Beta prior Beta(c,d) on rho
	double hpp_rho_c;
	double hpp_rho_d;

	double ibbs_loglik;
	double mc_loglik;
  
  // matrix of assignment probabilities for each locus to each state
  // used after MCMC is run for Stephens label switching algorithm
	arma::mat state_probs;

public:
  // constructors
	StickyHMM(int max_s,
           arma::uvec cts,
           arma::umat chr_begin_end,
           arma::vec m,
           double s);
  
  StickyHMM(int max_s,
            arma::uvec cts,
            arma::umat chr_begin_end,
            arma::vec m,
            double s,
            arma::uvec init_hidden_states,
            arma::vec init_state_rcn);
  
	StickyHMM(arma::umat chr_begin_end,
           arma::vec m,
           double s);
	
	// functions for simulating data
	void simulateData(int num_states,
                   int total_draws,
                   double self_trans_prob,
                   double gamma_a,
                   double gamma_scale);
	void simulateDataGivenStateCopyNum(int num_states,
                                    double self_trans_prob,
                                    arma::vec RCNs,
                                    int total_draws);
	void simulateTransitionMatrix(int num_states,
                               double self_trans_prob);
	void simulateHiddenStates();
	void simulateCopyNumbers(int num_states, 
                          double a,
                          double scale);
	void simulateCountsFromLociCopyNum(arma::vec lcns,
                                    int total_draws,
                                    arma::mat trans_mat,
                                    arma::rowvec init_dist,
                                    arma::uvec hs);
	void setLociCopyNumber();
	void calculateSimulatedDataLogLik();
	
	// getter functions
	arma::uvec getHiddenStates();
	arma::mat getTransMatrix();
	arma::rowvec getInitDist();
	arma::vec getLociCopyNumber();
	arma::vec getLociRelativeCopyNumber(const arma::vec &copy_number,
                                     const arma::vec &m);
	arma::vec getStateCopyNumber();
	arma::uvec getCounts();
	int max_s();
	double getIbbsll();
	double getMcll();

	// MCMC 
	void runMCMC(int iterations,
              arma::rowvec jump_size_prob,
              std::string sample_ref="",
              int thin=1,
              double gamma_a=3,
              double gamma_scale=1,
              double norm_sigma_s=0.15,
              double norm_sigma_l=0.5,
              double precision_sigma_s=50000,
              double precision_sigma_l=100000,
              double hpp_g_a=0.001,
              double hpp_g_b=0.001,
              double hpp_ak_a=0.01,
              double hpp_ak_b=0.0001,
              double hpp_r_c=1500,
              double hpp_r_d=1,
              double precision_p_shape = 3.874839,
              double precision_p_scale = 84594.12,
              double gg = -1, // gamma if passed as fixed value
              double rr = -1, // rho if passed as fixed value
              double apk = -1); // alpha + kappa if passed as fixed value
              
private:
  void countTransitions();
  void setTransMatGivenTransCounts();
  arma::mat computeMessages(int chr_index);
  arma::vec logLikLocusAllStates(const int &locus_index);
  void sampleStates(const int &chr_index, const arma::mat &messages);
  void sampleAuxiliaryVars();
  void sampleGlobalTransDist();
  arma::rowvec sampleFromDirichlet(const arma::rowvec &a);
  arma::rowvec sampleFromDirichlet_2(const arma::rowvec &a);
  void sampleTransitionMatrix();
  void GEM();
  double stateCopyNumberLogLik(const arma::uvec &loci_indices,
                               const double &CN);
  void sampleCopyNumber(const double &gamma_a,
                        const double &gamma_scale,
                        const double &norm_sigma_l,
                        const double &norm_sigma_s,
                        const arma::rowvec &jump_size_prob);
  void samplePrecision(const double &precision_sigma_s,
                       const double &precision_sigma_l,
                       const arma::rowvec &jump_size_prob);
  void sampleHyperParams();
  void calcLoglikCache();
  double stateOldCopyNumberLogLik(const arma::uvec &loci_indices,
                                  const int &state);
  void calculate_markov_chain_loglik();
  void calculate_IBBs_loglik();
};

#endif