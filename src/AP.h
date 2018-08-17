/**
 * Author: Sam Abujudeh
 * 
 * The Assignment Problem Class
 * 
 * Implemented as part of conliga for the purposes of
 * seeking the permutation of labels which minimuses the
 * cost function in the Stephens algorithm.
 * 
 * The Hungarian Algorithm (Munkres) is implemented here
 * and was based on description of the algorithm here
 * http://csclab.murraystate.edu/~bob.pilgrim/445/munkres.html
 * 
 * The Stephens algorithm is implemented in the LabelSwitch Class
 **/


#ifndef AP_H
#define AP_H

//#define ARMA_NO_DEBUG
#include <iostream>
#include <vector>
#include <armadillo>
#include <cfloat>

using namespace std;

class AP { //AP for Assignment Problem
	arma::mat cost;
	bool done;
	unsigned int step;
	unsigned int N;
	arma::umat M;
	arma::ivec rowCover;
	arma::ivec colCover;
	int cpath_0;
	int rpath_0;
	arma::imat path;

public:
	AP();
	AP(arma::mat c);
	arma::umat runMunkres();
	void set(arma::mat c);

private:
  // functions used by Munkres algorithm
	void step_one();
	void step_two();
	void step_three();
	void step_four();
	void step_five();
	void step_six();
	void find_noncovered_zero(int &row, int &col); //step four
	bool star_in_row(int &row); // step four
	void find_star_in_row(const int &row, int &col); // step four
	void find_star_in_col(const int &col, int &row); // step five
	void find_prime_in_row(const int &row, int &col); // step five
	void augment_path(const int &path_count); // step five
	void clear_covers(); 
	void erase_primes(); // step six
	void find_smallest(double &minval); // step six
};

#endif