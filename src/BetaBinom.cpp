#include "BetaBinom.h"     

using namespace Rcpp;

IndependentBBs::IndependentBBs() {}

IndependentBBs::IndependentBBs(arma::vec precs, arma::vec m)
/*
	Constructor used for simulating from a set of
	independent Beta Binomials (multiple samples)
*/
{
	means = m;
	precisions = precs;
	N_counts = means.size();
	N_samples = precisions.size();
	counts = arma::zeros<arma::umat>(N_counts, N_samples);
}

IndependentBBs::IndependentBBs(double prec, arma::vec m)
/*
	Constructor used for simulating from a set of
	independent Beta Binomials (one sample)
*/
{
	means = m;
	precisions.set_size(1);
	precisions[0] = prec;
	N_counts = means.size();
	N_samples = 1;
	counts = arma::zeros<arma::umat>(N_counts, N_samples);
}

IndependentBBs::IndependentBBs(arma::umat cnts)
/*
	Constructor used for fitting independent Beta Binomials
	given a matrix of counts for a set of samples
*/
{
	counts = cnts;
	N_counts = counts.n_rows;
	N_samples = counts.n_cols;
	precisions = arma::zeros<arma::vec>(N_samples);
	total_counts = arma::sum(counts,0); // sum cols
	calcLChooseCache();
	means = arma::zeros<arma::vec>(N_counts);
}

IndependentBBs::IndependentBBs(arma::uvec cnts, arma::vec m)
/*
	Constructor used for StickyHMM.
*/
{
	N_counts = cnts.size();
	N_samples = 1;
	counts.set_size(N_counts, N_samples);
	counts.col(0) = cnts;
	total_counts = arma::sum(counts,0); // sum cols
	calcLChooseCache();
	means = m;
	precisions = arma::zeros<arma::vec>(N_samples);
}

IndependentBBs::IndependentBBs(arma::uvec cnts, double prec, arma::vec m)
/*
	Constructor used for obtaining quantiles.
*/
{
	N_counts = cnts.size();
	N_samples = 1;
	counts.set_size(N_counts, N_samples);
	counts.col(0) = cnts;
	total_counts = arma::sum(counts,0); // sum cols
	//calcLChooseCache();
	means = m;
	precisions.set_size(1);
	precisions[0] = prec;
}


void IndependentBBs::calcLChooseCache()
{
	lchoose_cache.set_size(N_counts, N_samples);
	for (int i = 0; i < N_counts; ++i) {
		for (int j = 0; j < N_samples; ++j) {
			lchoose_cache(i, j) = Rf_lchoose(total_counts[j], counts(i, j));
		}
	}
}


void IndependentBBs::simulateCount(const int &sample, const int &locus, const int &N_draws)
{
	//double prec = precisionGivenMean(sample, means[locus]);
	//prec = 100000;
	
	double alpha = means[locus] * precisions[sample];
	double beta = (1 - means[locus]) * precisions[sample];
	double prob = Rf_rbeta(alpha, beta);
	counts(locus, sample) = Rf_rbinom(N_draws, prob);
}

void IndependentBBs::simulateCounts(arma::urowvec total_draws)
{
	total_counts = total_draws;
	for (int s = 0; s < N_samples; ++s) {
		for (int i = 0; i < N_counts; ++i) {
			simulateCount(s, i, total_draws[s]);
		}
	}
	total_counts = arma::sum(counts,0);
	calcLChooseCache(); // in order to calculate likelihood
}

arma::umat IndependentBBs::getCounts()
{
	return counts;
}

arma::uvec IndependentBBs::getSampleCounts(const int &sample)
{
	return counts.col(sample);
}

void IndependentBBs::calcMeansFromCopyNumber(const arma::vec &copy_number)
{
	arma::vec scaled_means = means % copy_number;
	double delta = arma::sum(means) / (arma::sum(scaled_means));
	means = delta * scaled_means;
}

arma::vec IndependentBBs::calcRelativeCopyNumFromCopyNumber(const arma::vec &copy_number,
															const arma::vec &m)
{
	means = m; // reset means because means are changed when simulating data
	arma::vec scaled_means = means % copy_number;
	double delta = arma::sum(means) / (arma::sum(scaled_means));
	return (delta * copy_number);
}

double IndependentBBs::dLogLik(const int &locus, const double &scaling_factor, const int &sample)
{
	double new_mean = means[locus]*scaling_factor;
	double alpha = calcAlpha(new_mean, sample);
	double beta = calcBeta(new_mean, sample);
	return (dLogBB(locus, sample, alpha, beta));
}

double IndependentBBs::dLogLikAllLoci(const arma::vec &scaling_factors, const double &precision, const int &sample)
/*
	If mean-var relationship is modelled instead:
	I will need to replace const double &precision with a vector of the
	coefs of the mean-var relationship.
*/
{
	/*
	Hack until I use a mean-var relationship
	I set the precision to the precision used to calculate the LogLik
	Then I set it back to what it was.
	*/
	double real_precision = getPrecision(means[0], sample);
	setPrecision(precision, sample);
	double ll = 0;
	for (int locus = 0; locus < N_counts; ++locus) {
		ll += dLogLik(locus, scaling_factors[locus], sample);
	}
	setPrecision(real_precision, sample);
	return ll;
}

double IndependentBBs::getPrecision(const double &mean, const int &sample)
/*
	I made it a function of the mean for future extension.
	However, the Precision is currently a constant for each sample.
	By default it outputs the precision of the first sample.
*/
{
	return precisions[sample];
}

void IndependentBBs::setPrecision(const double &precision, const int &sample)
{
	precisions[sample] = precision;
}

double IndependentBBs::calcAlpha(const double &mean, const int &sample)
{
	return (mean * getPrecision(mean, sample));
}

double IndependentBBs::calcBeta(const double &mean, const int &sample)
{
	double prec = getPrecision(mean, sample);
	return (prec - mean * prec);
}

double IndependentBBs::getAlpha(const int &locus, const int &sample)
{
/*
	If I cache the alphas I will need to change this to return
	the cached version instead.
*/
	return (calcAlpha(means[locus], sample));
}

double IndependentBBs::getMean(const int &locus)
{
	return means[locus];
}

arma::vec IndependentBBs::getMeans()
{
	return means;
}

void IndependentBBs::setMeans(arma::vec new_means)
{
	means = new_means;
}

void IndependentBBs::calcAlphasAndBetas()
{
	for (int s = 0; s < N_samples; ++s) {
		for (int i = 0; i < N_counts; ++i) {
			alphas(i, s) = calcAlpha(means[i], s);
			betas(i, s) = calcBeta(means[i], s);
		}
	}
}


double IndependentBBs::dLogBB(const int &locus, const int &sample,
							  const double &alpha, const double &beta)
{
	return (lchoose_cache(locus, sample) +
			Rf_lbeta(counts(locus, sample) + alpha, total_counts[sample] - counts(locus, sample) + beta) -
			Rf_lbeta(alpha, beta));
}

double IndependentBBs::dLogBBGivenMean(const int &locus, const int &sample,
									   const double &mean, const double &precision)
{
	//double prec = precisionGivenMean(sample, mean);
	double alpha = mean * precision;
	double beta = (1 - mean) * precision;
	return(dLogBB(locus, sample, alpha, beta));
}

double IndependentBBs::precisionLogLik(const int &sample, const double &precision)
/*
	If mean-var relationship is modelled instead:
	I will need to replace const double precision with a vector of the
	coefs of the mean-var relationship.
*/
{
	double ll = 0;
	for (int i = 0; i < N_counts; ++i) {
		ll += dLogBBGivenMean(i, sample, means[i], precision);
	}
	return ll;
}

void IndependentBBs::sampleMeans(const double &means_sigma)
{
	for (int i = 0; i < N_counts; ++i) {
		double new_mean = Rf_rnorm(means[i], means_sigma);
		if (new_mean >= 0) {
			double new_ll = 0;
			double old_ll = 0;
			for (int s = 0; s < N_samples; ++s) {
				new_ll += dLogBBGivenMean(i, s, new_mean, precisions[s]);
				old_ll += dLogBBGivenMean(i, s, means[i], precisions[s]);
			}
			double p = Rf_rexp(1);
			if (p > old_ll-new_ll) {
				means[i] = new_mean;
			}
		}
	}
}


void IndependentBBs::sampleMeans(const double &means_sigma, const arma::mat &mean_priors)
{
	for (int i = 0; i < N_counts; ++i) {
		double new_mean = Rf_rnorm(means[i], means_sigma);
		if (new_mean >= 0) {
			double new_ll = 0;
			double old_ll = 0;
			for (int s = 0; s < N_samples; ++s) {
				new_ll += dLogBBGivenMean(i, s, new_mean, precisions[s]);
				old_ll += dLogBBGivenMean(i, s, means[i], precisions[s]);

			}
			new_ll += Rf_dbeta(new_mean, mean_priors(0,i), mean_priors(1,i), 1);
			old_ll += Rf_dbeta(means[i], mean_priors(0,i), mean_priors(1,i), 1);

			double p = Rf_rexp(1);
			if (p > old_ll-new_ll) {
				means[i] = new_mean;
			}
		}
	}
}


/*
void IndependentBBs::samplePrecisionsCoefs(const arma::vec &precision_coefs_sigmas)
{
	for (int s = 0; s < N_samples; ++s) {
		arma::vec new_coefs(N_coefs);
		arma::vec old_coefs = precision_coefs.col(s);
		for (int i = 0; i < N_coefs; ++i) {
			new_coefs[i] = Rf_rnorm(precision_coefs[i], precision_coefs_sigmas[i]);
		}
		double new_ll = precisionLogLik(s, new_coefs);
		double old_ll = precisionLogLik(s, old_coefs);
		double p = Rf_rexp(1);
		if (p > old_ll-new_ll) {
			precision_coefs.col(s) = new_coefs;
		}
	}
}
*/

void IndependentBBs::samplePrecisions(const double &precision_sigma)
{
	for (int s = 0; s < N_samples; ++s) {
		double new_precision = Rf_rnorm(precisions[s], precision_sigma);
		if (new_precision > 0) {
			double new_ll = precisionLogLik(s, new_precision);
			double old_ll = precisionLogLik(s, precisions[s]);
			double p = Rf_rexp(1);
			if (p > old_ll-new_ll) {
				precisions[s] = new_precision;
			}
		} 
	}
}

void IndependentBBs::samplePrecisions(const double &precision_sigma,
	const arma::vec &precision_prior)
{
	for (int s = 0; s < N_samples; ++s) {
		double new_precision = Rf_rnorm(precisions[s], precision_sigma);
		if (new_precision > 0) {
			double new_ll = precisionLogLik(s, new_precision) +
				Rf_dgamma(new_precision, precision_prior[0], precision_prior[1], 1);
			double old_ll = precisionLogLik(s, precisions[s]) +
				Rf_dgamma(precisions[s], precision_prior[0], precision_prior[1], 1);
			double p = Rf_rexp(1);
			if (p > old_ll-new_ll) {
				precisions[s] = new_precision;
			}
		} 
	}
}

/*
double IndependentBBs::precisionGivenMean(const int &sample, const double &mean)
{
	double prec = precision_coefs(0, sample);
	double log_mean = log(mean);
	for (int i = 1; i < N_coefs; ++i) {
		prec += pow(log_mean, i) * precision_coefs(i, sample);
	}
	return exp(prec);
}
*/

// void IndependentBBs::inferParams(int iterations, double means_sigma, double precision_sigma, std::string run_ref)
// {
// 	// set means
// 	int all_counts = arma::accu(counts);
// 	for (int i = 0; i < N_counts; ++i) {
// 		for (int s = 0; s < N_samples; ++s) {
// 			means[i] += counts(i, s);
// 		}
// 		means[i] /= all_counts;
// 	}
// 
// 	ofstream means_out;
// 	ofstream precisions_out;
// 	means_out.open((run_ref + "_means.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
// 	precisions_out.open((run_ref + "_precisions.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
// 
// 	// need to put in a check that N_coefs is at least 2.
// 
// 	for (int s = 0; s < N_samples; ++s) {
// 		precisions[s] = 500000;
// 	}
// 	// randomly generate precision_coefs
// 
// 	for (int iter = 0; iter < iterations; ++iter) {
// 		Rcout << iter << std::endl;
// 		// update means
// 		sampleMeans(means_sigma);
// 		// update sample precisions
// 		samplePrecisions(precision_sigma);
// 
// 		//Rcout << means[0] << std::endl;
// 		for (int i = 0; i < N_counts; ++i) {
// 			means_out << std::setprecision(10) << means[i] << "\t";
// 		}
// 		means_out << "\n";
// 
// 		//means_out << std::fixed << std::setprecision(15) << means.t();
// 		for (int s = 0; s < N_samples; ++s) {
// 			precisions_out << precisions[s] << "\t";
// 		}
// 		precisions_out << "\n";
// 	}
// }


void IndependentBBs::inferParams(int iterations, double means_sigma, double precision_sigma,
		arma::mat mean_priors, arma::vec precision_prior,
		std::string run_ref)
{
	/* 	mean_priors is 2 by N_counts matrix
		precision_prior is a vector with [0]=shape and [1]=scale
	*/

	// set means
	int all_counts = arma::accu(counts);
	for (int i = 0; i < N_counts; ++i) {
		for (int s = 0; s < N_samples; ++s) {
			means[i] += counts(i, s);
		}
		means[i] /= all_counts;
	}

	ofstream means_out;
	ofstream precisions_out;
	means_out.open((run_ref + "_means.txt").c_str(), std::ofstream::out | std::ofstream::trunc);
	precisions_out.open((run_ref + "_precisions.txt").c_str(), std::ofstream::out | std::ofstream::trunc);

	// need to put in a check that N_coefs is at least 2.

	for (int s = 0; s < N_samples; ++s) {
		precisions[s] = Rf_rgamma(precision_prior[0], precision_prior[1]);
	}
	// randomly generate precision_coefs
	// randomly generate precision_coefs
	Rcout << "iteration: 1";
	for (int iter = 1; iter <= iterations; ++iter) {
	  if (iter % 1000 == 0) {
	    Rcout << iter;
	  } else {
	    if (iter % 100 == 0) {
	      Rcout << ".";
	    }
	  }

		// update means
		sampleMeans(means_sigma);
		// update sample precisions
		samplePrecisions(precision_sigma, precision_prior);

		//Rcout << means[0] << std::endl;
		for (int i = 0; i < N_counts; ++i) {
			means_out << std::setprecision(10) << means[i] << "\t";
		}
		means_out << "\n";

		//means_out << std::fixed << std::setprecision(15) << means.t();
		for (int s = 0; s < N_samples; ++s) {
			precisions_out << precisions[s] << "\t";
		}
		precisions_out << "\n";
	}
	Rcout << std::endl;
}

double IndependentBBs::dBetaBinom(double k, double N, double alpha, double beta)
{
  double log_dat = Rf_lchoose(N, k) + Rf_lbeta(k+alpha, N-k+beta) - 
    Rf_lbeta(alpha, beta);
  double prob = exp(log_dat);
  return prob;
}

double IndependentBBs::qBetaBinom(double p, double N, double alpha, double beta)
{
  double cummulate = 0;
  int indx_wanted;
  int i = 0;
  bool cont = true;
  // for(int i = 0; i < N; i++)
  while(cont)
  {
    cummulate += dBetaBinom(double(i), N, alpha, beta);
    if(cummulate > p){
      cont = false;
      indx_wanted = i;
    }
    i++;
  }
  return indx_wanted;
}

arma::mat IndependentBBs::getQuantiles(arma::vec q, int sample)
{
  //arma::mat quant_mat(alphas.size(), q.size());
  arma::vec alphas = means * precisions[sample];
  arma::mat quant_mat(alphas.size(), q.size());
  //double sum_alphas = sum(alphas);
  for (unsigned int i=0; i < alphas.size(); ++i) {
    for (unsigned int j=0; j < q.size(); ++j) {
      quant_mat(i,j) = qBetaBinom(q[j], total_counts[sample], alphas[i], precisions[sample]-alphas[i]);
    }
  }
  return quant_mat;
}

