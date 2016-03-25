/*
* randomVar.h
*/

#ifndef _RANDOMVAR_H_
#define _RANDOMVAR_H_
#include <time.h>

static int rnd_seed = -time(NULL);

namespace randomVar
{
	int setSeed(int seed);
	void set_rnd_seed (int s);
	double gammaln (const double z); /* returns ln(gamma(z)) */
	double factrl (const int n); /* returns n! */
	double logfactl (const int n); /* returns ln(n!) */
	double binomcoef (const int n, const int k); /* returns the binomial coefficient choose(n,k) */
	double binCoefTab (const int n, const int k); /* returns binomial coefficent choose(n,k) from table */
	double binomProb (const int n, const int k, const double p); /* Binomial(k;n,p) */
	double unif (); /* unif[0,1] */
	double wichmann_unif(int seed=0); /* U(0,1) */
	int poisson (double lamda); /* Pois(lamda)*/
	double fmax2(double x, double y); /* for rbeta */
	double fmin2(double x, double y); /* for rbeta */
	double rbeta(double aa, double bb); /* beta(alpha,beta)*/
	double uniform(int& seed = rnd_seed); /* a long period (> 2*10^18) generator for U[0,1] */
	double expo(double lamda); /* random deviate from exponential(lamda), lamda=rate parameter */
}

#endif /* _RANDOMVAR_H_ */
