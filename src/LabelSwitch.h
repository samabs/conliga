/**
 * Author: Sam Abujudeh
 * 
 * The Label Switch Class
 * 
 * Implemented as part of conliga for the purposes of
 * seeking the permutation of labels which minimuses the
 * cost function in the Stephens algorithm.
 * 
 * The Hungarian Algorithm (Munkres) is implemented in the
 * Assignment Problem (AP) class
 * and was based on description of the algorithm provided in:
 * http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
 * 
 * The Stephens algorithm is implemented in the LabelSwitch Class
 **/

#ifndef LABEL_SWITCH_H
#define LABEL_SWITCH_H

//#define ARMA_NO_DEBUG
#include <iostream>
#include <fstream>
#include <vector>
#include <armadillo>
#include <string>
#include "AP.h"
#include <typeinfo>


using namespace std;

class LabelSwitch { 
	arma::cube p; 
	arma::umat perm; 
	long long unsigned int K;
	long long unsigned int n;
	long long unsigned int m_iters;

public:
	LabelSwitch();

	LabelSwitch(std::string p_file,
             std::string n_states_file,
             long long unsigned int total_iterations, 
             long long unsigned int burn_in,
             long long unsigned int thin, 
             long long unsigned int K_states,
             long long unsigned int num_obs,
             bool map_max_state);
	
	arma::umat Stephens(int max_iter=100, double threshold=0.000001, bool identity=true);
};

#endif