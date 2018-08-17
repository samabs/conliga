#ifndef BETABINOM_H
#define BETABINOM_H


#include <RcppArmadillo.h>  
#include <iostream>
#include <vector>
#include "Sample.h"


using namespace std;


class IndependentBBs {
	int N_counts;
	int N_samples;
	arma::urowvec total_counts;
	arma::umat counts; // cols are samples, rows are loci
	arma::mat lchoose_cache; // cols are samples, rows are loci
	arma::vec precisions; // of size N_samples
	arma::mat alphas;
	arma::mat betas;
	arma::vec means; // of size N_counts

public:
  // constructors
	IndependentBBs();
	IndependentBBs(arma::uvec cnts, arma::vec m); // used by StickyHMM
	IndependentBBs(arma::vec precs, arma::vec m); // used for simulating data
	IndependentBBs(double prec, arma::vec m); // used for simulating data for one sample
	IndependentBBs(arma::umat cnts); 			  // used for fitting IBBs
	IndependentBBs(arma::uvec cnts, double prec, arma::vec m); // used for getting quantiles
	
	arma::umat getCounts();
	arma::uvec getSampleCounts(const int &sample=0);
	void simulateCount(const int &sample, const int &locus, const int &N_draws);
	void simulateCounts(arma::urowvec total_draws);
	double dLogBB(const int &locus, const int &sample,
				  const double &alpha, const double &beta);
	double dLogBBGivenMean(const int &locus, const int &sample,
					   	   const double &mean, const double &precision);
	double precisionLogLik(const int &sample, const double &precision);
	void sampleMeans(const double &means_sigma);

	void sampleMeans(const double &means_sigma, const arma::mat &mean_priors);

	void samplePrecisions(const double &precision_sigma);

	void samplePrecisions(const double &precision_sigma, const arma::vec &precision_prior);

	//void inferParams(int iterations, double means_sigma, double precision_sigma, std::string run_ref);

	void inferParams(int iterations, double means_sigma, double precision_sigma,
		arma::mat mean_priors, arma::vec precision_prior,
		std::string run_ref);

	void calcMeansFromCopyNumber(const arma::vec &copy_number);

	arma::vec calcRelativeCopyNumFromCopyNumber(const arma::vec &copy_number, const arma::vec &m);

	double getPrecision(const double &mean=0, const int &sample=0);
	void setPrecision(const double &precision, const int &sample=0);

	double calcAlpha(const double &mean, const int &sample=0);

	double calcBeta(const double &mean, const int &sample=0);

	double getMean(const int &locus);

	double getAlpha(const int &locus, const int &sample=0);

	arma::vec getMeans();
	
	void setMeans(arma::vec new_means);

	void calcAlphasAndBetas();

	double dLogLik(const int &locus, const double &scaling_factor=1.0, const int &sample=0);

	double dLogLikAllLoci(const arma::vec &scaling_factors, const double &precision, const int &sample=0);

	void calcLChooseCache();

	double dBetaBinom(double k, double N, double alpha, double beta);

	double qBetaBinom(double p, double N, double alpha, double beta);

	arma::mat getQuantiles(arma::vec q, int sample=0);

};

#endif

