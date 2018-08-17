#include "StickyHMM.h"
#include <google/profiler.h>


using namespace Rcpp;

// constructor used for HMM / copy number inference.
StickyHMM::StickyHMM(int max_s,
                     arma::uvec cts,
                     arma::umat chr_begin_end,
                     arma::vec m,
                     double s)
{
	max_states = max_s;
	chr_start_end = chr_begin_end;
	ibbs = IndependentBBs(cts, m);
	ibbs.setPrecision(s);
	num_chrs = chr_start_end.n_rows;
	num_loci = cts.size();
	num_loci_per_chr.set_size(num_loci);
	for (int i = 0; i < num_chrs; ++i) {
		num_loci_per_chr[i] = chr_start_end(i,1) - chr_start_end(i,0) + 1;
	}
	//hidden_states.set_size(max_states);
	//loci_RCNs.set_size(num_loci);
	hidden_states = arma::zeros<arma::uvec>(num_loci);
	state_RCNs = arma::zeros<arma::vec>(max_states);
	loci_RCNs = arma::zeros<arma::vec>(num_loci);
	trans_count = arma::zeros(max_states, max_states);
	trans_prob = arma::zeros(max_states, max_states);
	init_prob = arma::zeros<arma::rowvec>(max_states);
	w = arma::zeros<arma::vec>(max_states);
	m_bar_colsums = arma::zeros<arma::rowvec>(max_states);
	M = arma::zeros(max_states, max_states);
	M_bar = arma::zeros(max_states, max_states);
	loglik_cache = arma::zeros(max_states, num_loci);
	state_probs = arma::zeros(num_loci, max_states);
}

// constructor used for HMM / copy number inference.
StickyHMM::StickyHMM(int max_s,
                     arma::uvec cts,
                     arma::umat chr_begin_end,
                     arma::vec m,
                     double s,
                     arma::uvec init_hidden_states,
                     arma::vec init_state_rcn)
{
  max_states = max_s;
  chr_start_end = chr_begin_end;
  ibbs = IndependentBBs(cts, m);
  ibbs.setPrecision(s);
  num_chrs = chr_start_end.n_rows;
  num_loci = cts.size();
  num_loci_per_chr.set_size(num_loci);
  for (int i = 0; i < num_chrs; ++i) {
    num_loci_per_chr[i] = chr_start_end(i,1) - chr_start_end(i,0) + 1;
  }
  //hidden_states.set_size(max_states);
  hidden_states = arma::zeros<arma::uvec>(num_loci);
  hidden_states = init_hidden_states;
  state_RCNs = arma::zeros<arma::vec>(max_states);
  state_RCNs = init_state_rcn;
  loci_RCNs.set_size(num_loci);
  setLociCopyNumber();
  trans_count = arma::zeros(max_states, max_states);
  trans_prob = arma::zeros(max_states, max_states);
  init_prob = arma::zeros<arma::rowvec>(max_states);
  w = arma::zeros<arma::vec>(max_states);
  m_bar_colsums = arma::zeros<arma::rowvec>(max_states);
  M = arma::zeros(max_states, max_states);
  M_bar = arma::zeros(max_states, max_states);
  loglik_cache = arma::zeros(max_states, num_loci);
  state_probs = arma::zeros(num_loci, max_states);
}

// constructor used for simulating data.
// 
StickyHMM::StickyHMM(arma::umat chr_begin_end,
                     arma::vec m, double s)
{
	chr_start_end = chr_begin_end;
	ibbs = IndependentBBs(s, m);
	num_chrs = chr_start_end.n_rows;
	num_loci = chr_start_end(chr_start_end.n_rows-1, 1) + 1;
	arma::uvec num_loci_in_chr(num_loci);
	for (int i = 0; i < num_chrs; ++i) {
		num_loci_in_chr[i] = chr_start_end(i,1) - chr_start_end(i,0) + 1;
	}
	num_loci_per_chr = num_loci_in_chr;
	loci_RCNs.set_size(num_loci);
}


// simulate copy number profile and count data given:
// 1) number of true copy number states
// 2) number of total counts (total_draws)
// 3) the probability of self transition
// 4) the Gamma prior for over state RCN
// TODO: allow states and transition matrix to be simulated from Dirichlet Process
// TODO: allow for possibility of simulating inverse dispersion parameter from prior
void StickyHMM::simulateData(int num_states, int total_draws, double self_trans_prob,
	double gamma_a, double gamma_scale)
{
	simulateCopyNumbers(num_states, gamma_a, gamma_scale);
	simulateTransitionMatrix(num_states, self_trans_prob);
	simulateHiddenStates();
	setLociCopyNumber();
	simulateCountsFromLociCopyNum(loci_RCNs, total_draws, trans_prob, init_prob, hidden_states);
}

// the same as simulateData except, the state RCN is given
void StickyHMM::simulateDataGivenStateCopyNum(int num_states, double self_trans_prob,
	arma::vec RCNs, int total_draws)
/*
	RCNs: Note that it does not need to be relative copy number
 as this scaling is handled automatically
*/
{
	state_RCNs = RCNs;
	simulateTransitionMatrix(num_states, self_trans_prob);
	simulateHiddenStates();
	setLociCopyNumber();
	simulateCountsFromLociCopyNum(loci_RCNs, total_draws, trans_prob, init_prob, hidden_states);
}

// simulate RCN values for each state from the Gamma prior
// used in simulateData()
void StickyHMM::simulateCopyNumbers(int num_states, double a, double scale)
{
  state_RCNs.set_size(num_states);
  for (int i = 0; i < num_states; ++i) {
    state_RCNs[i] = Rf_rgamma(a, scale);
  }
}

// simulate transition matrix, given the self transition probability
// used in simulateData()
// TODO: allow simulation from Dirichlet Process
void StickyHMM::simulateTransitionMatrix(int num_states, double self_trans_prob)
{
  // set transition matrix.
  arma::mat real_trans_prob(num_states, num_states);
  double non_self_trans_prob = 1 - self_trans_prob;
  real_trans_prob.fill(non_self_trans_prob / (num_states-1));
  arma::vec self_trans(num_states);
  self_trans.fill(self_trans_prob);
  real_trans_prob.diag() = self_trans;
  trans_prob = real_trans_prob;
  
  init_prob = arma::ones<arma::rowvec>(num_states);
  init_prob = init_prob / num_states; // TODO: allow for non-uniform init_prob
}

// simulate hidden RCN states from Markov chain
// used in simulateData()
void StickyHMM::simulateHiddenStates()
{
  IntegerVector states = seq_len(trans_prob.n_rows)-1;
  hidden_states.set_size(num_loci);
  int loci_index = 0;
  for (int i = 0; i < num_chrs; ++i) {
    hidden_states[loci_index] = sample_int_fixed_prob(states); // add in init_prob
    loci_index += 1;
    for (int j = 1; j < num_loci_per_chr[i]; ++j) {
      arma::rowvec probs = trans_prob.row(hidden_states[loci_index-1]);
      hidden_states[loci_index] = sample_int(states, probs);
      loci_index += 1;
    }
  }
}

// given the loci RCN values, simulate the count data
// used in simulateData()
void StickyHMM::simulateCountsFromLociCopyNum(arma::vec lcns, int total_draws,
                                              arma::mat trans_mat, arma::rowvec init_dist,
                                              arma::uvec hs)
{
  trans_prob = trans_mat;
  init_prob = init_dist;
  hidden_states = hs;
  /*
   if we are simulating from loci copy number only,
   then we need the transition matrix, initial distribution
   which was used to generate the original hidden states in order to
   calculate the loglik. We also need to provide the hidden states.
   */

  // loci_RCNs needed for calculating loglik
  loci_RCNs = ibbs.calcRelativeCopyNumFromCopyNumber(lcns, ibbs.getMeans());
  arma::vec old_means = ibbs.getMeans(); // save means for controls
  ibbs.calcMeansFromCopyNumber(lcns); // this function changes the means of ibbs
  arma::urowvec N_draws(1);
  N_draws[0] = total_draws;
  ibbs.simulateCounts(N_draws);
  ibbs.setMeans(old_means); // set means back to WT in order to calculate likelihood function
  calculateSimulatedDataLogLik(); // to calculate ibbs_loglik and mc_loglik
  loci_RCNs = lcns; // set loci_RCN back to loci_CN so that we can get both CN and RCN
}

// used in simulateCountsFromLociCopyNum() to calculate the
// probability of the data given the parameters
void StickyHMM::calculateSimulatedDataLogLik() {
	Rcout << "mc" << std::endl;
	calculate_markov_chain_loglik();
	Rcout << "ibbs" << std::endl;
	calculate_IBBs_loglik();
}

// return the relative copy number, given the total copy number and the vector of means, m.
// used when simulating data
arma::vec StickyHMM::getLociRelativeCopyNumber(const arma::vec &copy_number, const arma::vec &m)
{
  // copy_number: the copy numbers for each locus
  return ibbs.calcRelativeCopyNumFromCopyNumber(copy_number, m);
}


// some getter functions
arma::uvec StickyHMM::getHiddenStates()
{
	return hidden_states;
}

arma::mat StickyHMM::getTransMatrix()
{
	return trans_prob;
}

arma::rowvec StickyHMM::getInitDist()
{
	return init_prob;
}

arma::vec StickyHMM::getLociCopyNumber()
{
	return loci_RCNs;
}

arma::vec StickyHMM::getStateCopyNumber()
{
  return state_RCNs;
}

arma::uvec StickyHMM::getCounts()
{
  return ibbs.getSampleCounts();
}

double StickyHMM::getMcll()
{
  return mc_loglik;
}

double StickyHMM::getIbbsll()
{
  return ibbs_loglik;
}

int StickyHMM::max_s()
{
  return max_states;
}

// for convenience, keep a vector of RCN for each locus
// rather than looking up the state RCN each time
void StickyHMM::setLociCopyNumber()
{
  for (int i = 0; i < num_loci; ++i) {
    loci_RCNs[i] = state_RCNs[hidden_states[i]];
  }
}

// used to calculate and store the output from the
// likelihood function to avoid calculating the likelihood
// function unnecessarily as it is computationally expensive
void StickyHMM::calcLoglikCache()
{
  for (int i = 0; i < num_loci; ++i) {
    loglik_cache.col(i) = logLikLocusAllStates(i);
  }
}

// used by calcLogLikCache() to calculate the log likelihood function
// for every locus in every state
arma::vec StickyHMM::logLikLocusAllStates(const int &locus_index)
{
  arma::vec logliks(max_states);
  double inf = std::numeric_limits<double>::infinity();
  for (unsigned int i=0; i < max_states; ++i) {
  	logliks[i] = ibbs.dLogLik(locus_index, state_RCNs[i]);
  }
  
  return(logliks);
}

// implementation of the stick breaking construction
// of the Dirichlet Process:
// sets the value of the variable beta when initialising MCMC
void StickyHMM::GEM() {
  arma::rowvec stick(max_states);
  double prev = 1;
  int i = 0;
  while (i < max_states) {
    stick[i] = Rf_rbeta(1, gamma) * prev;
    prev -= stick[i];
    i++;
  }
  beta = stick;
}

// compute messages for each Markov chain
// this is usually for each chomosome arm
// TODO: this can produce NaNs if chain starts in
//        unlikely place - need a fix.
//        Current solution is just to get user to 
//        run the MCMC again manually.
//        Could allow MCMC to restart a defined number of times
//        or start the chain from a more likely place
arma::mat StickyHMM::computeMessages(int chr_index)
{
	int loci_index;
	unsigned int nrows = num_loci_per_chr[chr_index];
 	unsigned int ncols = max_states;
  	arma::mat messages = arma::zeros(nrows, ncols);
  	messages.row(nrows-1).fill(1);
  	for (int r = nrows-1; r > 0; --r) {
  		loci_index = chr_start_end(chr_index, 0) + r;
  		messages.row(r-1) = (trans_prob * (exp(loglik_cache.col(loci_index)) % messages.row(r).t())).t();
  		messages.row(r-1) /= messages.row(r-1).max();
  		
  		//if (messages.row(r-1).has_nan()) {
  			// TODO: need to handle this:
  			// sometimes, if the MCMC starts in a very unlikely place,
  			// runMCMC() can fail.
  			// we could either: 1) start chain again randomly
  			//                  2) start chain off in more likely place
  			// Rcout << "locus: " << loci_index << std::endl;
  			// Rcout << state_RCNs.t() << std::endl;
  			// Rcout << loglik_cache.col(loci_index).t() << std::endl;
  			// Rcout << messages.row(r-1) << std::endl;
  			// Rcout << messages.row(r) << std::endl;
  			// Rcout << loglik_cache.col(loci_index+1).t() << std::endl;
  			// Rcout << messages.row(r+1) << std::endl;
  		//}

  	}
  	if (messages.has_nan()) {
  	  Rcout << "Warning! NaNs detected, MCMC will fail due to starting in an unlikely place..." << std::endl;
  	}
  	return(messages);
}


// used when initialising MCMC to set the
// initial transition matrix.
void StickyHMM::setTransMatGivenTransCounts()
{
  // set transition matrix (pi)
	trans_prob = arma::zeros(max_states, max_states);
	for (int i = 0; i < max_states; ++i) {
		if (sum(trans_count.row(i))==0) {
		  //double self_trans = alpha * beta[i] + kappa / (alpha + kappa);
		  //trans_prob.row(i).fill((1-self_trans)/max_states);
		  //trans_prob(i,i) = self_trans;
			trans_prob(i,i) = 1;
		} else {
			trans_prob.row(i) = trans_count.row(i) / sum(trans_count.row(i));
		}
		
	}
  
  // set initial distribution (pi^0)
	init_prob.zeros();
	arma::uvec first_loci_indices = chr_start_end.col(0);
 	arma::uvec first_loci_hs = hidden_states.elem(first_loci_indices);
  	arma::rowvec first_loci_hs_counts(max_states);
  	first_loci_hs_counts.zeros();
  	for (int i = 0; i < max_states; ++i) {
  		arma::uvec loci_indices_in_state = arma::find(first_loci_hs == i);
  		first_loci_hs_counts[i] = loci_indices_in_state.size();
  	}
  	init_prob = first_loci_hs_counts / sum(first_loci_hs_counts);
}

// given the current hidden state assignments,
// count the state transitions in each chromosome arm
void StickyHMM::countTransitions() {
	trans_count = arma::zeros(max_states, max_states);
	for (int i = 0; i < num_chrs; ++i) {
		for (int j = chr_start_end(i,0) + 1; j <= chr_start_end(i,1); ++j) {
			trans_count(hidden_states[j-1], hidden_states[j]) += 1;
		}
	}
}

// Sample the hidden states for a chromosome arm given the messages that were
// calculated by computeMessages() function
// GIBBS
void StickyHMM::sampleStates(const int &chr_index, const arma::mat &messages)
{

	IntegerVector states = seq_len(max_states)-1;
  // get first locus in chromosome arm so we access the correct loci
	int chr_offset = chr_start_end(chr_index, 0);
	
	// sample state for first locus in chromosome arm
	arma::rowvec logprobs(max_states);
	for (int j = 0; j < max_states; ++j) {
		logprobs[j] = log(init_prob[j]) + log(messages(0, j)) +
		loglik_cache(j, chr_offset);
	}

	arma::rowvec probs = arma::exp(logprobs-arma::max(logprobs));
	probs = probs / arma::sum(probs);
	state_probs.row(chr_offset) = probs;
	hidden_states[chr_offset] = sample_int(states, probs);

  // calculate probability of assignment to each state for each locus
	for (int l = 1; l < num_loci_per_chr[chr_index]; ++l) {
		arma::rowvec logprobs(max_states);
		unsigned int k = hidden_states[l+chr_offset];
		for (int j = 0; j < max_states; ++j) {
			logprobs[j] = log(messages(l, j)) +
						  log(trans_prob(hidden_states[l+chr_offset-1], j)) +
						  loglik_cache(j, l+chr_offset);

		}
		
		// convert logprobs to probs
		arma::rowvec probs = arma::exp(logprobs-arma::max(logprobs));
		probs = probs / arma::sum(probs);

		// save probability of assignments to states for later output
		// this is for resolving label switching using Stephens algorithm
		state_probs.row(l+chr_offset) = probs;
		
		// sample the states using the calculated probabilities
		hidden_states[l+chr_offset] = sample_int(states, probs);
	}
}


// sample auxiliary variables: M, w and M bar
// used in runMCMC()
// these auxiliary variables are required to sample the hyperparameters
void StickyHMM::sampleAuxiliaryVars()
{
    // sample M
  	M = arma::zeros(max_states, max_states);
  	arma::mat P = arma::zeros(max_states, max_states);
  	P.each_row() += beta;
  	P = P * alpha;
  	P.diag() += kappa;
  	
  	for (int j = 0; j < max_states; ++j) {
  	  for (int k = 0; k < max_states; ++k) {
  	    for (int n = 0; n < trans_count(j,k); ++n) {
  	      if (Rf_rbinom(1, P(j,k) / (n + P(j,k)))) {
  	        M(j,k) += 1;
  	      }
  	    }
  	  }
  	}
  	
  	// ensure rho is set
  	rho = kappa / (alpha + kappa);
  	
  	// sample w
  	w.set_size(beta.size());
  	for (unsigned int i = 0; i < beta.size(); ++i) {
    	w[i] = Rf_rbinom(M(i, i), rho / (rho + beta[i] * (1 - rho)));
  	}
  	
  	// set M_bar - used for sampling beta
  	M_bar = M;
  	M_bar.diag() -= w;
  	m_bar_colsums = arma::sum(M_bar, 0);
}


// sample beta in MCMC
void StickyHMM::sampleGlobalTransDist()
{
	arma::rowvec a = (gamma / max_states) + m_bar_colsums;
	beta = sampleFromDirichlet(a);
}


arma::rowvec StickyHMM::sampleFromDirichlet(const arma::rowvec &a)
{
	arma::rowvec y(a.size());
	for (int i = 0; i < a.size(); ++i) {
		y[i] = Rf_rgamma(a[i], 1);
	}
	y = y / arma::sum(y);
	if (y.has_nan()) {
		y = sampleFromDirichlet_2(a);
	}
	return y;
}

// Called by sampleFromDirichlet:
// Used when values of y are very small causing overflow (NaNs) when dividing by the sum of y
arma::rowvec StickyHMM::sampleFromDirichlet_2(const arma::rowvec &a)
{
	arma::rowvec y(a.size());
	y.zeros();
	y[0] = Rf_rbeta(a[0], accu(a.subvec(1,a.size()-1)));
	double phi = 0;
	for (int i = 0; i < a.size()-1; ++i) {
		phi = Rf_rbeta(a[i], accu(a.subvec(i+1,a.size()-1)));
		y[i] = (1-accu(y)) * phi;
	}
	y[a.size()-1] = 1-accu(y);
	return y;
}

// This function samples:
// 1) the initial distribution
// 2) the transition matrix
void StickyHMM::sampleTransitionMatrix()
{
  // Transition matrix
  arma::mat new_pi = arma::zeros(max_states, max_states);
  new_pi.each_row() = alpha * beta;
  new_pi += trans_count;
  new_pi.diag() += kappa;
  
  // sample transition matrix
  for (int i = 0; i < max_states; ++i) {
    new_pi.row(i) = sampleFromDirichlet(new_pi.row(i));
  }
  trans_prob = new_pi;
  
  // initial distribution (pi^0)
  // sample initial distribution using counts from first loci
  arma::uvec first_loci_indices = chr_start_end.col(0);
  arma::uvec first_loci_hs = hidden_states.elem(first_loci_indices);
  arma::urowvec first_loci_hs_counts(max_states);
  first_loci_hs_counts.zeros();
  for (int i = 0; i < max_states; ++i) {
    arma::uvec loci_indices_in_state = arma::find(first_loci_hs == i);
    first_loci_hs_counts[i] = loci_indices_in_state.size();
  }
  
  // sample initial distribution
  init_prob = sampleFromDirichlet((alpha * beta) + first_loci_hs_counts);
}

// used in Metropolis-Hastings step to update state RCN
// given the new loci assigned to the state,
// calculate the likelihood function using the proposed value of RCN
double StickyHMM::stateCopyNumberLogLik(const arma::uvec &loci_indices,
                                        const double &CN)
{
  double loglik = 0;
  for (int i = 0; i < loci_indices.size(); ++i) {
    int locus = loci_indices[i];
    loglik += ibbs.dLogLik(locus, CN);
  }
  return loglik;
}

// used in Metropolis-Hastings step to update state RCN
// given the new loci assigned to the state,
// calculate the likelihood function using the previous value of RCN
// Note we can use the stored cache values...
double StickyHMM::stateOldCopyNumberLogLik(const arma::uvec &loci_indices,
                                           const int &state)
{
  // create uvec to use submat function.
  arma::uvec ss(1);
  ss[0] = state;
  double loglik = arma::accu(loglik_cache.submat(ss, loci_indices));
  return loglik;
}

// Metropolis-Hastings step to update state RCN
void StickyHMM::sampleCopyNumber(const double &gamma_a, const double &gamma_scale,
                                 const double &norm_sigma_l, const double &norm_sigma_s,
                                 const arma::rowvec &jump_size_prob)
{
  for (int k = 0; k < max_states; ++k) {
    // select loci assigned to state k
    arma::uvec loci_indices = arma::find(hidden_states == k);
    double new_CN;
    if (loci_indices.size() > 0) { // if loci assigned to state...
      IntegerVector jump_choice = seq_len(2)-1;
      // choose to make a large or small jump
      int large_jump = sample_int(jump_choice, jump_size_prob);
      if (large_jump) {
        new_CN = state_RCNs[k] + Rf_rnorm(0, norm_sigma_l);
      }
      else {
        new_CN = state_RCNs[k] + Rf_rnorm(0, norm_sigma_s);
      }
      if (new_CN >= 0) { // proposed RCN must be 0 or greater (otherwise reject proposal)
        double loglik_old = stateOldCopyNumberLogLik(loci_indices, k) +
          Rf_dgamma(state_RCNs[k], gamma_a, gamma_scale, 1);
        double loglik_new = stateCopyNumberLogLik(loci_indices, new_CN) +
          Rf_dgamma(new_CN, gamma_a, gamma_scale, 1);
        double p = Rf_rexp(1);
        if (p > loglik_old-loglik_new) {
          state_RCNs[k] = new_CN; // accept proposal
        }
      }
    }
    else { // if no loci in state, sample from gamma prior
      state_RCNs[k] = Rf_rgamma(gamma_a, gamma_scale);
    }
  }
  // update loci_RCNs to reflect new state_RCNs
  setLociCopyNumber();	
}

// sample inverse overdispersion parameter (Note we call this precision here)
void StickyHMM::samplePrecision(const double &precision_sigma_s,
								const double &precision_sigma_l,
								const arma::rowvec &jump_size_prob)
{
	double loglik = 0;
	IntegerVector jump_choice = seq_len(2)-1;
	double proposed_precision;
	// choose large or small jump
	int large_jump = sample_int(jump_choice, jump_size_prob);
	
	// propose new value
	if (large_jump) {
		proposed_precision = Rf_rnorm(ibbs.getPrecision(), precision_sigma_l);
	}
	else {
		proposed_precision = Rf_rnorm(ibbs.getPrecision(), precision_sigma_s);
	}

	if (proposed_precision > 0) { // proposal must be greater than 0, otherwise reject

	  
		double p = Rf_rexp(1);
	  /*
	   * TODO: could save these values in the cache so don't have to calculate again
	   * when calculating the cache at the beginning of an iteration.
	   * This would only save us having to calculate 1/max_s of the values so
	   * might not speed things up that much
	   * 
	   */
		double old_loglik = ibbs.dLogLikAllLoci(loci_RCNs, ibbs.getPrecision()) +
			Rf_dgamma(ibbs.getPrecision(), precision_prior_shape, precision_prior_scale, 1);
		double new_loglik = ibbs.dLogLikAllLoci(loci_RCNs, proposed_precision) +
			Rf_dgamma(proposed_precision, precision_prior_shape, precision_prior_scale, 1);

		if (p > old_loglik-new_loglik) {
			ibbs.setPrecision(proposed_precision); // accept proposal
			ibbs_loglik = new_loglik; // save value for output
		} else { // otherwise reject proposal
			ibbs_loglik = old_loglik;
		}
	}
}

// sample the hyperparameters:
// gamma, alpha and kappa
void StickyHMM::sampleHyperParams() {
	arma::vec trans_count_rowsums = arma::sum(trans_count, 1);
	double sum_M = arma::accu(M);
	double sum_M_bar = arma::accu(M_bar);
	arma::rowvec rr(max_states);
	rr.ones();
	arma::rowvec ss(max_states);
	ss.zeros();
	double temp_shape;
	double temp_rate;

	// ################################## (alpha+kappa) ##################################
	double alpha_plus_kappa = alpha + kappa;
	
	if (!apk_fixed) {
	  for (int j = 0; j < max_states; ++j) {
	    if (trans_count_rowsums[j] > 0) {
	      ss[j] = Rf_rbinom(1, trans_count_rowsums[j] / (trans_count_rowsums[j] + alpha + kappa));
	      rr[j] = Rf_rbeta(alpha + kappa + 1, trans_count_rowsums[j]);
	    }
	  }
	  
	  temp_shape = hpp_aplusk_a + sum_M - arma::sum(ss);
	  temp_rate = hpp_aplusk_b - arma::sum(arma::log(rr));
	  
	  //Rf_rgamma parameterised by shape and scale, where scale = 1/rate
	  alpha_plus_kappa = Rf_rgamma(temp_shape, 1/temp_rate);
	  
	}

	// ################################## (gamma) ##################################
	
	if (!gamma_fixed) {
	  arma::urowvec temp_gt = m_bar_colsums > 0;
	  double K_bar = arma::accu(temp_gt);
	  
	  // sample zeta
	  int zeta = Rf_rbinom(1, sum_M_bar / (sum_M_bar + gamma));
	  
	  // sample eta
	  double eta = Rf_rbeta(gamma + 1, sum_M_bar);
	  
	  // sample gamma
	  temp_shape = hpp_gamma_a + K_bar - zeta;
	  temp_rate = hpp_gamma_b - log(eta);
	  gamma = Rf_rgamma(temp_shape, 1/temp_rate); // Gamma(shape, scale)
	}

	// ################################## (rho) ##################################

	if (!rho_fixed) {
	  double sum_w = arma::accu(w);
	  rho = Rf_rbeta(sum_w + hpp_rho_c, sum_M - sum_w + hpp_rho_d);
	}

	// ############################# (update kappa and alpha) ##############################
	
	kappa = rho * alpha_plus_kappa;
	alpha = alpha_plus_kappa * (1 - rho);
	
}


// calculate probability of hidden states
void StickyHMM::calculate_markov_chain_loglik() {
	mc_loglik = 0;
	arma::mat log_trans_prob = log(trans_prob);
	arma::rowvec log_init_prob = log(init_prob);
	for (int c = 0; c < num_chrs; ++c) {
		unsigned int chr_start = chr_start_end(c, 0);
		unsigned int chr_end = chr_start_end(c, 1);
		mc_loglik += log_init_prob[hidden_states[chr_start]];
		for (unsigned int locus = chr_start; locus < chr_end; ++locus) {
			mc_loglik += log_trans_prob(hidden_states[locus], hidden_states[locus+1]);
		}
	}
}

void StickyHMM::calculate_IBBs_loglik() {
	ibbs_loglik = ibbs.dLogLikAllLoci(loci_RCNs, ibbs.getPrecision());
}

// Main algorithm - initialise and run MCMC
void StickyHMM::runMCMC(int iterations,
                        arma::rowvec jump_size_prob,
                        std::string sample_ref,
                        int thin,
                        double gamma_a, // prior on RCN
                        double gamma_scale, // prior on RCN
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
                        double apk) // alpha + kappa if passed as fixed value

{
	/*
	Random start

	Inferring all:
		hidden states
		RCN
		precision
		hyperparameters
	*/
	
	// ############### INITIALISE MCMC ###############

	// set hyper-param prior parameters
	hpp_gamma_a = hpp_g_a;
	hpp_gamma_b = hpp_g_b;
	hpp_aplusk_a = hpp_ak_a;
	hpp_aplusk_b = hpp_ak_b;
	hpp_rho_c = hpp_r_c;
	hpp_rho_d = hpp_r_d;

	//set gamma prior for precision
	precision_prior_shape = precision_p_shape;
	precision_prior_scale = precision_p_scale;

	// output MCMC parameters
	ofstream params_out;
	params_out.open((sample_ref + "_params.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	params_out << "iterations\tthin\tmax_states\tnum_loci\tnum_chr_arms\tjump_size_prob_small\tjump_size_prob_large";
	params_out << "\tRCN_gamma_a\tRCN_gamma_scale\tRCN_norm_sigma_s\tRCN_norm_sigma_l\tprecision_norm_sigma_s\tprecision_norm_sigma_l";
	params_out << "\thpp_gamma_a\thpp_gamma_b\thpp_aplusk_a\thpp_aplusk_b\thpp_rho_c\thpp_rho_d\tprecision_p_shape\tprecision_p_scale" << std::endl;
	params_out << iterations << "\t" << thin << "\t" << max_states << "\t" << num_loci << "\t" << num_chrs << "\t"<< jump_size_prob[0] << "\t" << jump_size_prob[1];
	params_out << "\t" << gamma_a << "\t" << gamma_scale << "\t" << norm_sigma_s << "\t" << norm_sigma_l << "\t" << precision_sigma_s << "\t" << precision_sigma_l;
	params_out << "\t" << hpp_gamma_a << "\t" << hpp_gamma_b << "\t" << hpp_aplusk_a << "\t" << hpp_aplusk_b << "\t" << hpp_rho_c << "\t" << hpp_rho_d << "\t" << precision_prior_shape << "\t" << precision_prior_scale;
	params_out.close();
	params_out.clear();

	gamma_fixed = false;
	rho_fixed = false;
	apk_fixed = false;
	
	if (gg != -1) { // gamma give as fixed
	  gamma_fixed = true;
	  gamma = gg;
	} else { // sample gamma
	  gamma = Rf_rgamma(hpp_gamma_a, 1/hpp_gamma_b); // b = inverse scale
	}
	
	if (rr != -1) { // rho given as fixed
	  rho_fixed = true;
	  rho = rr;
	} else { // sample rho
	  rho = Rf_rbeta(hpp_rho_c, hpp_rho_d);
	}
	
	double alpha_plus_kappa;
	if (apk != -1) { // alpha_plus_kappa given as fixed
	  apk_fixed = true;
	  alpha_plus_kappa = apk;
	} else { // sample alpha_plus_kappa
	  alpha_plus_kappa = Rf_rgamma(hpp_aplusk_a, 1/hpp_aplusk_b);
	}
	
	// set values of kappa and alpha
	kappa = rho * alpha_plus_kappa;
	alpha = alpha_plus_kappa - kappa;
	
	// draw precision from prior
	ibbs.setPrecision(Rf_rgamma(precision_prior_shape, precision_prior_scale));

	//define output files and clear them
	ofstream hs_out;
	ofstream cn_out;
	ofstream p_out;
	ofstream tm_out;
	ofstream alpha_out;
	ofstream gamma_out;
	ofstream kappa_out;
	ofstream mc_loglik_out;
	ofstream ibbs_loglik_out;
	ofstream beta_out;
	ofstream state_probs_out;
	ofstream n_pop_states_out;
	hs_out.open((sample_ref + "_hs.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	cn_out.open((sample_ref + "_cn.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	p_out.open((sample_ref + "_p.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	alpha_out.open((sample_ref + "_alpha.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	gamma_out.open((sample_ref + "_gamma.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	kappa_out.open((sample_ref + "_kappa.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	mc_loglik_out.open((sample_ref + "_mcll.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	ibbs_loglik_out.open((sample_ref + "_ibbsll.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	beta_out.open((sample_ref + "_beta.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	state_probs_out.open((sample_ref + "_probs.bin").c_str(), ios::out | ios::binary);
	n_pop_states_out.open((sample_ref + "_nstates.bin").c_str(), ios::out | ios::binary);
	
	// for outputting state_probs to binary file:
	long long unsigned int n_bytes = sizeof(state_probs(0,0));
	long long unsigned int n_loci = num_loci;
	long long unsigned int n_states = max_states;

	// vector of state indices
	IntegerVector states = seq_len(max_states)-1;
	
	// if not provided initial hidden state assignments
	// then randomly draw them
	if (arma::sum(hidden_states) == 0) { 
	  // randomly assign loci to states
	  hidden_states = sample_with_replacement_fixed_prob(states, num_loci);
	  
	  // draw the copy number from a gamma dist.
	  simulateCopyNumbers(max_states, gamma_a, gamma_scale);
	  
	  // set loci copy num to state copy num
	  setLociCopyNumber();
	}

	// stick breaking process
	// returns rowvec beta
	GEM();

	// set trans_count given hidden states
	countTransitions();

	// set M
	M = arma::zeros(max_states, max_states);

	setTransMatGivenTransCounts();

	// calc loglik of IBBs and markov chain
	calculate_IBBs_loglik();
	calculate_markov_chain_loglik();

	// output initial values of chain
	for (int locus_index=0; locus_index < num_loci-1; ++locus_index)
	{
		//output hidden states
		hs_out << hidden_states[locus_index] << "\t";
	}
	// output last locus...
	hs_out << hidden_states[num_loci-1] << "\n";

	for (int state_index=0; state_index < max_states-1; ++state_index)
	{
		//output state RCN
		cn_out << state_RCNs[state_index] << "\t";
		//output beta
		beta_out << beta[state_index] << "\t";
	}
	// output last state
	cn_out << state_RCNs[max_states-1] << "\n";
	beta_out << beta[max_states-1] << "\n";

	// output precision
	p_out << ibbs.getPrecision() << "\n";

	// output hyperparameters
	alpha_out << alpha << "\n";
	gamma_out << gamma << "\n";
	kappa_out << kappa << "\n";

	// output logliks
	mc_loglik_out << mc_loglik << "\n";
	ibbs_loglik_out << ibbs_loglik << "\n";

	//output transition matrix
	// ostringstream ss;
	// ss << 0;
	// tm_out.open((sample_ref + "_tm" + ss.str() + ".txt").c_str());
	// tm_out << trans_prob;
	// tm_out.close();
	// tm_out.clear();

	
	// ############### START ITERATIONS ###############
	
	// 1) for each chromosome:
	//   i) compute messages
	//   ii) sample states
	// 2) sample auxiliary variables
	// 3) update global transition distribution
	// 4) sample transition matrix
	// 5) sample "emission": hidden state copy number
	// 6) sample precision of Beta Binomials
	// 7) sample hyperparameters (optionally)
	
	Rcout << "starting MCMC for sample: " << sample_ref << std::endl;

	for (int iter = 1; iter <= iterations; ++iter) {

		calcLoglikCache();

		for (int c = 0; c < num_chrs; ++c) {
			// step 1.i)
			arma::mat messages = computeMessages(c);

			// step 1.ii)
			sampleStates(c, messages);
		}
		// update trans_count after sampling states
		// TODO: update transition when sampling each state
		countTransitions();

		// step 2)
		sampleAuxiliaryVars();
		
		// step 3)
		sampleGlobalTransDist();

		// step 4)
		sampleTransitionMatrix();

		// step 5)
		sampleCopyNumber(gamma_a, gamma_scale, norm_sigma_s, norm_sigma_l, jump_size_prob);

		// step 6) 
		samplePrecision(precision_sigma_s, precision_sigma_l, jump_size_prob);
		
		// ibbs_loglik set by samplePrecision
		calculate_markov_chain_loglik();

		// step 7)
		sampleHyperParams();
		//gamma = 1;

		// output to files when  
		if (iter % thin == 0) {
			for (int locus_index=0; locus_index < num_loci-1; ++locus_index)
			{
				//output hidden states
				hs_out << hidden_states[locus_index] << "\t";
			}
			// output last locus...
			hs_out << hidden_states[num_loci-1] << "\n";

			for (int state_index=0; state_index < max_states-1; ++state_index)
			{
				//output state RCN
				cn_out << state_RCNs[state_index] << "\t";
				// output beta
				beta_out << beta[state_index] << "\t";
			}
			// output last state
			cn_out << state_RCNs[max_states-1] << "\n";
			beta_out << beta[max_states-1] << "\n";

			// output precision
			p_out << ibbs.getPrecision() << "\n";

			// output hyperparameters
			alpha_out << alpha << "\n";
			gamma_out << gamma << "\n";
			kappa_out << kappa << "\n";

			// output logliks
			mc_loglik_out << mc_loglik << "\n";
			ibbs_loglik_out << ibbs_loglik << "\n";

			//output transition matrix
			// ostringstream ss;
			// ss << iter;
			// tm_out.open((sample_ref + "_tm" + ss.str() + ".txt").c_str());
			// tm_out << trans_prob;
			// tm_out.close();
			// tm_out.clear();

			arma::uvec pop_states = arma::unique(hidden_states); // in ascending order
			arma::mat pop_state_probs = state_probs.cols(pop_states);
			unsigned int n_pop_states = pop_states.size();
			// output loci state probabilities
			state_probs_out.write((char*) pop_state_probs.memptr(), n_bytes * n_pop_states * n_loci );
			// output number of populated states in the interation
			n_pop_states_out.write(reinterpret_cast<const char*>(&n_pop_states), sizeof(unsigned int));

		}

	}
}


