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

#include "LabelSwitch.h"

using namespace std;

LabelSwitch::LabelSwitch() {}

LabelSwitch::LabelSwitch(std::string p_file,
                         std::string n_states_file,
                         long long unsigned int total_iterations, 
                         long long unsigned int burn_in,
                         long long unsigned int thin, 
                         long long unsigned int K_states,
                         long long unsigned int num_obs,
                         bool map_num_state) {
  // map_num_state is true if using MAP method and false if using MAX method
  
  // this constructor creates a LabelSwitch object and reads in the assignment probabilities
  
	/*use long long unsigned int to prevent overflow
	for the number of bytes in large binary files. */
	
	n = num_obs;
	m_iters = (total_iterations - burn_in) / thin;
	long long unsigned int n_bytes = sizeof(double);
	
	ifstream in_n_file(n_states_file.c_str(), ios::in | ios::binary | ios::ate);
	
	cout << "reading in state assignment probability data..." << std::endl;

	in_n_file.seekg(0);
	vector<unsigned int> n_pop_states;
	n_pop_states.resize(total_iterations / thin);
	in_n_file.read(reinterpret_cast<char *>(&n_pop_states[0]), (total_iterations / thin) * sizeof(n_pop_states[0]));
	arma::uvec n_states(&n_pop_states[0], total_iterations / thin);

	ifstream in_p_file(p_file.c_str(), ios::in | ios::binary | ios::ate);
	streampos f_size = in_p_file.tellg();
	long long unsigned int sum_burned_ps = 0;
	if (burn_in > 0) {
		cout << "burn in > 0..." << std::endl;
		cout << "number of iterations to burn: " << (burn_in / thin) << std::endl;
		sum_burned_ps = arma::accu(n_states.subvec(0, (burn_in/thin)-1)) * n;
		cout << "number of probabilities to burn: " << sum_burned_ps << std::endl;
	}

	n_states = n_states.subvec(burn_in / thin, (total_iterations / thin)-1);
	cout << "n_states after burn in: " << n_states.size() << std::endl;
	unsigned int map_state;
	unsigned int map_count = 0;
	if (map_num_state) {
		arma::uvec unique_n_states = arma::unique(n_states);
		
		for (unsigned int i = 0; i < unique_n_states.size(); ++i) {
			arma::uvec positions = arma::find(n_states == unique_n_states[i]);
			unsigned int state_count = positions.size();
			if (state_count > map_count) {
				map_count = state_count;
				map_state = unique_n_states[i];
			}
		}
		K = map_state;
		p = arma::zeros<arma::cube>(n, K, map_count);
		cout << "MAP number of states: " << map_state << std::endl;
		cout << "Number of iterations to be kept: " << map_count << std::endl;
	} else {
		K = arma::max(n_states);
		p = arma::zeros<arma::cube>(n, K, m_iters);
	}

	long long unsigned int start_pos = n_bytes * sum_burned_ps;
    in_p_file.seekg(start_pos);

	if (map_num_state) {
		cout << "reading in probabilities only from " << map_count << " iterations with the MAP number of states..." << std::endl;
		unsigned int p_idx = 0;
		for (unsigned int m = 0; m < m_iters; ++m) {
			vector<double> buffer(K*n, 0.0);
			if (n_states[m] == map_state) {
				//cout << p_idx << " of " << map_count << std::endl;
				in_p_file.read(reinterpret_cast<char *>(&buffer[0]), K * n * sizeof(buffer[0]));
				p.slice(p_idx) = arma::mat(&buffer[0], n, K, false);
				p_idx++;
			} else {
				in_p_file.ignore(n_states[m] * n * sizeof(buffer[0]));
			}
		}
		m_iters = p_idx;
		cout << m_iters << std::endl;
		perm.set_size(m_iters, K);
	} else {
		for (unsigned int m = 0; m < m_iters; ++m) {
			long long unsigned int temp = n_states[m];
			vector<double> buffer(temp*n, 0.0);
			in_p_file.read(reinterpret_cast<char *>(&buffer[0]), temp * n * sizeof(buffer[0]));
			arma::mat temp_mat = arma::mat(&buffer[0], n, temp, false);
			p.slice(m).submat(0, 0, n-1, temp-1) = temp_mat; 
		}
		perm.set_size(m_iters, K);
	}

	arma::urowvec states = arma::linspace<arma::urowvec>(0, (K-1), K);
	perm.each_row() = states;
	cout << "finished reading in data." << std::endl;
}


arma::umat LabelSwitch::Stephens(int max_iter, double threshold, bool identity) {

  /* TODO: at the moment identity doesn't have any use
   * the idea would be to allow the algoritm to start 
   * not in the identity condition.  This would be useful
   * if the state switching solution didn't work and we could
   * try running it again but not from the identity condition.
   */
	cout << perm.row(0) << std::endl;
	cout << "replace 0s and 1s in p" << std::endl;

	p.elem(arma::find(p < 0.000001)).fill(0.000001);

	p.elem(arma::find(p > (1-0.000001))).fill(1-0.000001);

	cout << "make rows of p equal to 1." << std::endl;
	arma::vec p_sum(K);
	for (unsigned int m = 0; m < m_iters; ++m) {
		p_sum = sum(p.slice(m), 1);
		for (unsigned i = 0; i < n; ++i) {
			p.slice(m).row(i) /= p_sum[i];
		}
	}

	cout << "starting Stephens algorithm" << std::endl;
	int iter = 0;
	double criterion = 99;
	AP ap;
	arma::umat ap_res(K, K);
	arma::mat temp = arma::zeros(n, K);
	arma::mat q = arma::zeros(n, K);
	arma::mat cost = arma::zeros(K, K);
	arma::mat log_q = arma::zeros(n, K);
	arma::mat log_p_m = arma::zeros(n, K);
	arma::uvec index(1);
	double previous = -99;
	double current = 0;
	
	while ((criterion > threshold) && (iter < max_iter)) {
		iter++;
		cout << "iteration: " << iter;
		cout << ", calculating q";
		q = arma::zeros(n, K);
		
		for (unsigned int m = 0; m < m_iters; ++m) {
		  q += p.slice(m).cols(perm.row(m));
		}
		
		q = q / m_iters;
		cout << ", making cost matrices";
		log_q = log(q);
		
		for (unsigned int m = 0; m < m_iters; ++m) {
			log_p_m = log(p.slice(m));
		  
			for (unsigned int j=0; j < K; ++j) {
				temp = log_p_m;
				temp.each_col() -= log_q.col(j);
				temp %= p.slice(m);
				cost.row(j) = sum(temp, 0);
			}
			
			ap.set(cost);
			ap_res = ap.runMunkres();
			
			for (unsigned int c = 0; c < K; ++c) {
				index = arma::find(ap_res.col(c) == 1, 1);
				perm(m, c) = index[0];
			}
			
			perm.row(m) = arma::sort_index(perm.row(m)).t();
			current += accu(cost);
		}
		
		criterion = abs(previous - current);
		cout << ", criterion: " << criterion << std::endl;
		previous = current;
		current = 0;
	}
	
	if (criterion > threshold) {
		cout << "Max iterations exceeded" << std::endl;
	}
	
	return(perm);
}