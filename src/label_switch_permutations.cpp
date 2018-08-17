/**
 * Make the Label Switching functions available for use by conliga
 **/

#ifndef LS_ATTR_CPP
#define LS_ATTR_CPP

#include "LabelSwitch.h"
#include "AP.h"

using namespace std;


// [[Rcpp::export]]
arma::umat get_state_permutations(std::string p_file, std::string n_states_file,
                                  int total_iterations, int burn_in, int thin, int K_states, int num_obs,
                                  bool map_num_state, int stephens_max_iters)
{
  // Function called when processing the MCMC to perform the Stephens label switching algorithm
  
  LabelSwitch ls(p_file,
                 n_states_file,
                 total_iterations,
                 burn_in,
                 thin,
                 K_states,
                 num_obs,
                 map_num_state);
  
  arma::umat perm = ls.Stephens(stephens_max_iters, 0.000001, true);
  return (perm);
}

#endif