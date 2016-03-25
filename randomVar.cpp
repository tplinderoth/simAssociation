/*
* randomVar.cpp
*/

#include "randomVar.h"
#include "Matrix.h"
#include <math.h>
#include <vector>
#include <cstdlib>
#include <cstdio> // rbeta
#include <cfloat> // rbeta

#define expmax	(DBL_MAX_EXP * M_LN2)/* = log(DBL_MAX) */

// FUNCTION DEFINITIONS

int randomVar::setSeed(int seed)
{
	int static z_rndu;
	z_rndu = 170*(seed%178) + 137;
	return z_rndu;
}

double randomVar::gammaln (const double z)
{
/*
* Numerical Recipes 2nd ed. pg 219
* returns ln(gamma(z)) for z > 0
*/

	if (z <= 0)
	{
		fprintf(stderr, "invalid value %f given to randomVar::gammaln\n", z);
		return -1;
	}

	int j;
	double x, y, tmp, ser;
	static const double cof[6]={76.18009172947146, -86.50532032941677,
		24.01409824083091, -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5};

	y=x=z;
	tmp=x+5.5;
	tmp -= (x+0.5)*log(tmp);
	ser=1.000000000190015;
	for (j=0; j<6; ++j)
		ser += cof[j]/++y;
	return -tmp+log(2.5066282746310005*ser/x);
}

double randomVar::factrl (const int n)
{
/*
* Numerical Recipes 2nd ed. pg 219
* returns n!
*/

	static int ntop=4;
	const static int overflow_thresh = 33;
	static double a[overflow_thresh]={1.0, 1.0, 2.0, 6.0, 24.0};
	int j;

	if (n < 0)
	{
		fprintf(stderr, "Negative value %i used in randomVar::factrl\n", n);
		return -1;
	}
	if (n > overflow_thresh - 1)
		return exp(gammaln(n+1.0));
	while (ntop<n)
	{
		j=ntop++;
		a[ntop]=a[j]*ntop;
	}
	return a[n];
}

double randomVar::logfactl (const int n)
{
/*
* Numerical Recipes 2nd ed. page 220
* returns ln(n!)
*/

	const static int tablesz=101;
	static double a[tablesz];

	if (n < 0)
	{
		fprintf(stderr, "Negative value %i used in randomVar::logfactrl\n", n);
		return -1;
	}
	else if (n <= 1.0)
		return 0.0;
	else if (n < tablesz)
		return (a[n] != 0.0 ? a[n] : (a[n]=gammaln(n+1.0)));
	else
		return gammaln(n+1.0);
}

double randomVar::binomcoef (const int n, const int k)
{
/*
* Numerical Recipes 2nd ed. page 220
* return the binomial coefficient choose(n,k)
*/

	if (n < 0.0 || k < 0.0)
	{
		fprintf(stderr, "Negative value used for randomVar::binomcoef: choose(%i,%i)\n",n,k);
		return -1;
	}
	else if (n < 1.0 || k > n)
		return 0.0;
	else
		return floor(0.5+exp(logfactl(n)-logfactl(k)-logfactl(n-k)));

}

double randomVar::binCoefTab (const int n, const int k)
{
/*
* returns the binomial coefficient choose(n,k)
* use when many different binomial coefficients need to be computed
* increased speed at the expense of memory usage for storing a table
*/

	static Matrix<double> coef(101, 101);

	if (n < 1.0 || k > n)
		return 0.0;

	static unsigned int nidx = static_cast<unsigned int> (n);
	static unsigned int kidx = static_cast<unsigned int> (k);

	if (nidx < coef.rown() && kidx < coef.coln())
		return(coef[nidx][kidx] != 0.0 ? coef[nidx][kidx] : (coef[nidx][kidx]=binomcoef(n,k)));
	else
		return binomcoef(n,k);
}

double randomVar::binomProb (const int n, const int k, const double p)
{
	if ( (k < 0.0 || k > n ) || p > 1.0 || p < 0.0)
	{
		fprintf(stderr, "invalid arguments to randomVar::binomProb n=%i, k=%i, p=%f\n",n,k,p);
		return -1;
	}
	else if (p == 0.0)
	{
		if (k > 0)
			return 0.0;
		else
			return 1.0;
	}
	else if (p == 1.0)
		if (k < n)
			return 0.0;
		else
			return 1.0;
	else
		return exp(logfactl(n) - logfactl(k) -logfactl(n-k) + static_cast<double>(k)*log(p) + static_cast<double>(n-k)*log(1.0-p));
}

double randomVar::unif ()
{
	return (static_cast<double>(rand()) / (RAND_MAX) );
}

double randomVar::wichmann_unif(int seed)
{
	/*U(0,1): AS 183: Appl. Stat. 31:188-190
	Wichmann BA & Hill ID.  1982.  An efficient and portable pseudo-random number generator.  Appl. Stat. 31:188-190 x, y, z are any numbers in the range 1-30000.  Integer operation up to 30323 required.
	Suggested to me by Z. Yang who also provided me with the source code used here.
	*/
	static int x_rndu=11, y_rndu=23;
	static int z_rndu;
	if (seed!=0)
		z_rndu = setSeed(seed);
	else
		z_rndu=137;
	double r;
	x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
	y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
	z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
	if (x_rndu<0) x_rndu+=30269;
	if (y_rndu<0) y_rndu+=30307;
	if (z_rndu<0) z_rndu+=30323;
	r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
	return (r-(int)r);
}

// poisson draws random number from Poisson(lambda)
int randomVar::poisson (double lamda)
{
	double u = unif();
	double p = exp(-lamda);
	int i = 0;
	int x = i;

	while (u > p)
	{
		p += (lamda * p) / (i + 1);
		i++;
		x = i;
	}

	return x;
}

double randomVar::rbeta(double aa, double bb)
{
	/*
	Taken from R project 2.13.2
	Thorfinn@binf.ku.dk oct 07 2011
	This should properly be implemented better using gsl or boost.
	we should cite the R project ....
	 */

	/* Reference:
	 * R. C. H. Cheng (1978).
	 * Generating beta variates with nonintegral shape parameters.
	 * Communications of the ACM 21, 317-322.
	 * (Algorithms BB and BC)
	 */

    double a, b, alpha;
    double r, s, t, u1, u2, v, w, y, z;

    int qsame;
    /* FIXME:  Keep Globals (properly) for threading */
    /* Uses these GLOBALS to save time when many rv's are generated : */
    static double beta, gamma, delta, k1, k2;
    static double olda = -1.0;
    static double oldb = -1.0;

    if (aa <= 0. || bb <= 0. || (isinf(aa) && isinf(bb)))
      fprintf(stderr,"Warning: alpha and beta parameters passed to randomVar::rbeta are %f and %f\n",aa,bb);

    if (isinf(aa))
    	return 1.0;

    if (isinf(bb))
    	return 0.0;


    /* Test if we need new "initializing" */
    qsame = (olda == aa) && (oldb == bb);
    if (!qsame) { olda = aa; oldb = bb; }

    a = fmin2(aa, bb);
    b = fmax2(aa, bb); /* a <= b */
    alpha = a + b;

#define v_w_from__u1_bet(AA) 			\
	    v = beta * log(u1 / (1.0 - u1));	\
	    if (v <= expmax) {			\
		w = AA * exp(v);		\
		if(isinf(w)) w = DBL_MAX;	\
	    } else				\
		w = DBL_MAX


    if (a <= 1.0) {	/* --- Algorithm BC --- */

	/* changed notation, now also a <= b (was reversed) */

	if (!qsame) { /* initialize */
	    beta = 1.0 / a;
	    delta = 1.0 + b - a;
	    k1 = delta * (0.0138889 + 0.0416667 * a) / (b * beta - 0.777778);
	    k2 = 0.25 + (0.5 + 0.25 / delta) * a;
	}
	/* FIXME: "do { } while()", but not trivially because of "continue"s:*/
	for(;;) {
	    u1 = unif();
	    u2 = unif();
	    if (u1 < 0.5) {
		y = u1 * u2;
		z = u1 * y;
		if (0.25 * u2 + z - y >= k1)
		    continue;
	    } else {
		z = u1 * u1 * u2;
		if (z <= 0.25) {
	 	    v_w_from__u1_bet(b);
		    break;
		}
		if (z >= k2)
		    continue;
	    }

	    v_w_from__u1_bet(b);

	    if (alpha * (log(alpha / (a + w)) + v) - 1.3862944 >= log(z))
		break;
	}
	return (aa == a) ? a / (a + w) : w / (a + w);

    }
    else {		/* Algorithm BB */

	if (!qsame) { /* initialize */
	    beta = sqrt((alpha - 2.0) / (2.0 * a * b - alpha));
	    gamma = a + 1.0 / beta;
	}
	do {
	    u1 = unif();
	    u2 = unif();

	    v_w_from__u1_bet(a);

	    z = u1 * u1 * u2;
	    r = gamma * v - 1.3862944;
	    s = a + r - w;
	    if (s + 2.609438 >= 5.0 * z)
		break;
	    t = log(z);
	    if (s > t)
		break;
	}
	while (r + alpha * log(alpha / (b + w)) < t);

	return (aa != a) ? b / (b + w) : w / (b + w);
    }
}

double randomVar::fmax2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? y : x;
}

double randomVar::fmin2(double x, double y)
{
#ifdef IEEE_754
	if (ISNAN(x) || ISNAN(y))
		return x + y;
#endif
	return (x < y) ? x : y;
}


double randomVar::expo(double lamda)
{
	return -log(uniform())/lamda;
}


double randomVar::uniform(int &seed)
{
	/*
	 * "ran2" from Numerical Recipes 2nd ed. page 286
	 * Returns a uniform random deviate from [0,1]
	 * Call uniform() with the initializing 'seed' < 0; thereafter, do not alter 'seed' between successive deviates in a sequence.
	 */

	const int IM1=2147483563, IM2=2147483399;
	const int IA1=40014, IA2=40692, IQ1=53668, IQ2=52774;
	const int IR1=12211, IR2=3791, NTAB=32, IMM1=IM1-1;
	const int NDIV=1+IMM1/NTAB;
	const double EPS=3.0e-16, RNMX=1.0-EPS, AM=1.0/(double)IM1;
	static int seed2=123456789, iy=0;
	static std::vector<int> iv(NTAB);
	int j,k;
	double temp;

	if (seed <= 0) // initialize
	{
		seed=(seed==0 ? 1 : -seed);
		seed2=seed;
		for(j=NTAB+7;j>=0;--j)
		{
			k=seed/IQ1;
			seed=IA1*(seed-k*IQ1)-k*IR1;
			if (seed < 0) seed += IM1;
			if (j < NTAB) iv[j] = seed;
		}
		iy=iv[0];
	}
	k=seed/IQ1;
	seed=IA1*(seed-k*IQ1)-k*IR1;
	if (seed < 0) seed += IM1;
	k=seed2/IQ2;
	seed2=IA2*(seed2-k*IQ2)-k*IR2;
	if (seed2 < 0) seed2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-seed2;
	iv[j]=seed;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

void randomVar::set_rnd_seed (int s)
{
	rnd_seed=s;
}
