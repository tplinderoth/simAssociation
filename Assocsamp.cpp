/*
 * Assocsamp.cpp
 */

#include "Assocsamp.h"
#include "randomVar.h"
#include <cstdio>
#include <cstdlib>
#include <time.h>
#include <iostream>
#include <math.h>

void Pileup::allocate_read_storage (int depth)
{
	reads.resize(depth);
	qual.resize(depth);
}

Assocsamp::Assocsamp ()
	: _ntreatments(2),
	  _npools(0),
	  _counts(NULL),
	  _minorid(-1),
	  _minorf(NULL),
	  _fst(0.0),
	  _fail(0)
{
	srand((unsigned)time(NULL));
}

Assocsamp::~Assocsamp ()
{
	int i,j;

	// clean up vector of pool data
	if (!_pools.empty())
	{
		for (i=0; i<_npools; ++i)
		{
			for (j=0; j<_ntreatments; ++j)
				delete [] _pools[i][j];
			delete [] _pools[i];
		}
		_npools = 0;
	}

	// release minor allele frequency info storage
	if (_minorf)
		delete [] _minorf;
	_minorf=NULL;

	// release _counts memory
	if (_counts)
	{
		for (i=0; i<_ntreatments; ++i)
			delete [] _counts[i];
		delete [] _counts[i];
		_counts=NULL;
	}

	_ntreatments = 0;
}

void Assocsamp::initialize_pool (int npool, int ntreatment)
{
	int i,j,k;

	if (npool < 1)
	{
		fprintf(stderr,"Invalid number of pools %i passed to Assocsamp::initalize_pool\n", npool);
		_fail = 1;
		return;
	}
	if (ntreatment < 1)
	{
		fprintf(stderr,"Invalid number of treatments %i passed to Assocsamp::initalize_pool\n", ntreatment);
		_fail = 1;
		return;
	}

	_pools.resize(npool);

	// reserve memory and initialize
	for (i=0; i<npool; ++i)
	{
		_pools[i] = new int* [_ntreatments];
		for (j=0; j<_ntreatments; ++j)
		{
			_pools[i][j] = new int [2];
			for (k=0; k<2; ++k)
				_pools[i][j][k] = 0;
		}
	}

	_npools = npool;
	_ntreatments = ntreatment;

	initMinorfreq();
	initAlleleCounts();
}

void Assocsamp::initMinorfreq ()
{
	_minorf = new double[_ntreatments+1];
	for (int i=0; i<=_ntreatments; ++i)
		_minorf[i] = 0.0;
}

void Assocsamp::initAlleleCounts ()
{
	int i,j;

	_counts = new int* [_ntreatments];

	for (i=0; i<_ntreatments; ++i)
	{
		_counts[i] = new int [2];
		for (j=0; j<2; ++j)
			_counts[i][j] = 0;
	}
}

int& Assocsamp::allele1_count (int pool, int treatmentid)
{
	if (treatmentid > _npools - 1)
	{
		fprintf(stderr,"treatment ID %i out of bounds in call to Assocsamp::allele1_count\n",treatmentid);
		exit(EXIT_FAILURE);
	}
	return _pools[pool][treatmentid][0];
}

int& Assocsamp::allele2_count (int pool, int treatmentid)
{
	if (treatmentid > _npools - 1)
	{
		fprintf(stderr,"treatment ID %i out of bounds in call to Assocsamp::allele2_count\n",treatmentid);
		exit(EXIT_FAILURE);
	}
	return _pools[pool][treatmentid][1];
}


int Assocsamp::allele1_count (int pool, int treatmentid) const
{
	if (treatmentid > _npools - 1)
	{
		fprintf(stderr,"treatment ID %i out of bounds in call to Assocsamp::allele1_count\n",treatmentid);
		exit(EXIT_FAILURE);
	}
	return _pools[pool][treatmentid][0];
}

int Assocsamp::allele2_count (int pool, int treatmentid) const
{
	if (treatmentid > _npools - 1)
	{
		fprintf(stderr,"treatment ID %i out of bounds in call to Assocsamp::allele2_count\n",treatmentid);
		exit(EXIT_FAILURE);
	}
	return _pools[pool][treatmentid][1];
}

int Assocsamp::quantitative_pool (int nind, int npools, double f [], int treatmentid, int ntreatments, int nsubpop, double fst)
{
	/*
	 * nind: diploid pool size
	 * npools: number of pools
	 * f: frequency of causative allele; f=P(allele0) in each population, has length=nsubpop
	 * treatmentid: index of treatment that the values of 'f' belong to (i.e 0,1,2...,etc.)
	 * ntreatments: number of total treatments in analysis
	 * nsubpop: number of populations to sample from
	 * fst: Fst between the sampled populations
	 *
	 */

	static int genoc [3]; // [00, 01, 11]
	int n = nind*npools;
	int remain = n;
	double draw;
	int i,j;

	if (_pools.empty())
	{
		initialize_pool(npools, ntreatments);
		if (_fail)
			return -1;
	}
	if (_pools.size() != static_cast<unsigned int>(npools))
	{
		_fail=1;
		fprintf(stderr,"Number of pools passed to Assocsamp::quantitative_pool differs from _pools container size\n");
		return -1;
	}

	if (nsubpop > 1)
	{
		_fst = fst;
		if (fst == 0.0)
			fprintf(stderr,"warning: sampling %i populations with an Fst of zero\n",nsubpop);
		sampleMultipops(n, fst, nsubpop, f, genoc, 0); // sample genotypes from multiple populations
	}
	else
	{
		if (fst > 0.0)
			fprintf(stderr,"warning: nonzero Fst has no effect for one population\n");
		genGenotypes(f[0], n, genoc); // sample genotypes from 1 population
	}

	// log allele count information
	_counts[treatmentid][0] = 2*genoc[0]+genoc[1];
	_counts[treatmentid][1] = 2*genoc[2]+genoc[1];

	// pool individuals
	for (i=0; i<npools; ++i) // figure out how to make sampling for pools look more random
	{
		for (j=0; j<nind; ++j)
		{
			draw = randomVar::uniform();
			if (draw <= genoc[0]/remain)
			{
				allele1_count(i,treatmentid) += 2;
				--genoc[0];
			}
			else if (draw >= (1-genoc[2]/remain))
			{
				allele2_count(i,treatmentid) += 2;
				--genoc[2];
			}
			else
			{
				++allele1_count(i,treatmentid);
				++allele2_count(i,treatmentid);
				--genoc[1];
			}
			--remain;
		}
	}

	return 0;
}

void Assocsamp::sampleMultipops (int n, double fst, int npops, double f [], int totalgeno [], bool isancestral)
{
	/*
	 * n: diploid sample size
	 * fst: fst among the populations
	 * npops: number of subpopulations to sample from
	 * f: allele frequency in each population
	 * totalgeno: array to hold the counts of genotypes for the total population sample
	 * isancestral: are the supplied allele frequencies for the ancestral population?
	 */

	static int subgeno[3]; // [00, 01, 11]
	int subn[npops];
	int i,j;

	for (i=0; i<3; ++i)
		totalgeno[i]=0;

	// sample equal number of individuals among the subpopulations
	for (i=0; i<npops; ++i)
		subn[i]=n/npops;
		j=0;
		int sz = arraysum(subn);
		while (sz != n)
		{
			if (sz > n)
				--subn[j];
			else
				++subn[j];
			j = j+1 >= npops ? 0 : j+1;
			sz=arraysum(subn);
		}

		for (i=0; i<npops; ++i)
		{
			// draw allele frequencies from Balding-Nichols Distribution
			if (isancestral)
				f[i]=randomVar::rbeta((((1.0-fst)/fst)*f[0]), (((1.0-fst)/fst)*(1.0-f[0])));
			// sample genotypes
			genGenotypes(f[i], subn[i], subgeno);
			for (j=0; j<3; ++j)
				totalgeno[j] += subgeno[j];
		}
}

int Assocsamp::arraysum(int a [])
{
	int sum=0;
	for (unsigned int i=0; i<sizeof(a)/sizeof(a[0]); ++i)
		sum += a[i];
	return sum;
}

void Assocsamp::genGenotypes (double f, int n, int geno [])
{
	/*
	 * f = frequency of '0' allele
	 * n: diploid sample size
	 * geno: [00, 01, 11]
	 */

	static double expected [3];
	static int count [3];
	int sum = 0;
	int i,j;

	expected[0]=f*f*n;
	expected[1]=2*f*(1-f)*n;
	expected[2]=(1-f)*(1-f)*n;

	for (i=0; i<3; ++i)
	{
		count[i] = expected[i];
		geno[i] = (expected[i] - count[i]) < 0.5 ? count[i] : count[i]+1;
		sum += geno[i];
	}

	while (sum != n)
	{
		if (randomVar::uniform() < 0.5) // avoid systematic bias
		{
			i=0;
			j=2;
		}
		else
		{
			i=2;
			j=0;
		}
		if ( (expected[i]-count[i]) == 0.5 )
			geno[i] = sum > n ? geno[i]-1 : geno[i]+1;
		else if ( expected[j]-count[j] == 0.5)
			geno[j] = sum > n ? geno[j]-1 : geno[j]+1;
		else
			geno[1] = sum > n ? geno[1]-1 : geno[1]+1;

		sum = 0;
		for (i=0; i<3; ++i)
			sum += geno[i];
	}
}

void Assocsamp::setMinorfreq ()
{
	int n0, n1, hapsz;
	double a0f, a1f;
	int i;

	n0=n1=0;
	for (i=0; i<_ntreatments; ++i)
	{
		n0 += _counts[i][0];
		n1 += _counts[i][1];
	}
	hapsz = n0+n1;
	a0f=static_cast<double>(n0)/hapsz;
	a1f=static_cast<double>(n1)/hapsz;
	_minorid = a0f < a1f ? 0 : 1;
	_minorf[_ntreatments] = a0f < a1f ? a0f : a1f;

	for (i=0; i<_ntreatments; ++i)
	{
		_minorf[i] = static_cast<double>(_counts[i][_minorid])/(_counts[i][0]+_counts[i][1]);
	}
}

void Assocsamp::treatGenos (double f, int n, double dominance, bool iscase, int geno [])
{
	/*
	 * f: frequency of causative allele 'a' in population
	 * n: number of sampled diploid individuals
	 * dominance: degree of dominance of the causative allele
	 * iscase: sampling case (vs control) individuals?
	 * geno: [number aa, number aA, number AA]
	 */

	double paa = f*f;
	double paA = 2*f*(1-f);
	double pAA = (1-f)*(1-f);
	double ppheno = iscase ? paa+paA : pAA+paA; // figure out how to use dominance to model quantitative traits
	static double expected [3];
	static int count[3];
	int sum;
	int i;

	if (iscase)
	{
		expected[0] = dominance == 1.0 ? n*(paa/ppheno) : 1.0;
		expected[1] = dominance == 1.0 ? n*(paA/ppheno) : 0.0;
		expected[2] = 0.0;
	}
	else
	{
		expected[0] = 0.0;
		expected[1] = dominance == 1.0 ? 0.0 : n*(pAA/ppheno);
		expected[2] = dominance == 1.0 ? 1.0 : n*(paA/ppheno);
	}

	for(i=0; i<3; ++i)
	{
		count[i] = expected[i];
		geno[i] = (expected[i]-count[i]) < 0.5 ? count[i] : count[i]+1;
	}

	// make sure genotype counts add up to the sample size
	sum=0;
	i = iscase ? 0 : 2;
	sum = geno[1]+geno[i];
	while( sum != n)
	{
		if (sum < n)
			++geno[i];
		else
			--geno[i];
		sum = geno[1]+geno[i];
	}
}

int Assocsamp::makePileup (double avgcov, double error, double maxqual, double offsetq, std::vector<Pileup>* pilevec)
{
	int i,j,k;
	const int extra = 20;
	char reference, alternate;
	unsigned int n = _npools*_ntreatments;
	std::vector<Pileup>* pile = pilevec ? pilevec : &_pile;

	if (_pools.empty())
	{
		fprintf(stderr,"No data to create pileup from in call to Assocsamp::makePileup\n");
		return -1;
	}

	setMinorfreq();
	if (pile->size() != n)
		pile->resize(n);

	setAlleleID(&reference, &alternate);
	k=0;
	for (i=0; i<_ntreatments; ++i)
	{
		for (j=0; j<_npools; ++j)
		{
			if ((*pile)[k].reads.empty())
				(*pile)[k].allocate_read_storage(avgcov+extra);
			(*pile)[k].ref = reference;
			(*pile)[k].alt = alternate;
			if (_pools[j][i][0]+_pools[j][i][1] < 2) // treatment has no data
				(*pile)[k].cov = 0;
			else
				(*pile)[k].cov = randomVar::poisson(avgcov);
			pileReads(minorfreq(i), &(*pile)[k]);
			pileQuality(error, maxqual, offsetq, &(*pile)[k]);
			++k;
		}
	}

	return 0;
}

void Assocsamp::setAlleleID (char* ref, char* alt)
{
	// set reference allele
	*ref=drawBase();

	// set alternate allele
	*alt=*ref;
	while (*alt == *ref)
		*alt=drawBase();
}

char Assocsamp::drawBase ()
{
	double draw = randomVar::unif();
	char b;

	if (draw <= 0.25)
		b='A';
	else if (draw > 0.25 && draw <= 0.50)
		b='C';
	else if (draw > 0.5 && draw <= 0.75)
		b='G';
	else
		b='T';

	return b;
}

void Assocsamp::pileReads (double minorfreq, Pileup* pile)
{
	char r;
	double allele_draw, strand_draw;

	if (pile->cov == 0) // missing data
		pile->reads[0] = '*';
	else
	{
		for (int i=0; i<pile->cov; ++i)
		{
			allele_draw = randomVar::uniform();
			strand_draw = randomVar::uniform();
			if (allele_draw <= minorfreq)
				r = strand_draw < 0.5 ? tolower(pile->alt) : pile->alt;
			else
				r = strand_draw < 0.5 ? ',' : '.';

			if (i < static_cast<int>(pile->reads.size()))
				pile->reads[i] = r;
			else
				pile->reads.push_back(r);
		}
	}
}

void Assocsamp::pileQuality (double rate, double maxq, double offset, Pileup* pile)
{
	int qscore;
	char qascii;

	if (pile->cov == 0) // missing data
		pile->qual[0] = '*';
	else
	{
		for (int i=0; i<pile->cov; ++i)
		{
			qscore = -10.0*log10(randomVar::expo(1.0/rate)); // generate error from an exponential distribution
			if (qscore > maxq) qscore = maxq; // make sure quality score does not exceed maximum
			qascii = qscore+offset;
			if (i < static_cast<int>(pile->qual.size()))
				pile->qual[i] = qascii;
			else
				pile->qual.push_back(qascii);
			doError(qscore, pile->ref, pile->alt, pile->reads[i]);
		}
	}
}

void Assocsamp::doError (int qscore, char ref, char alt, char& read)
{
	double err = pow(10.0, static_cast<double>(-qscore)/10.0);
	double draw = randomVar::uniform();
	int isforward;
	char r, errorbase;

	if (draw <= err)
	{
		// determine strand
		r = toupper(read);
		if (read == '.')
		{
			isforward = 1;
			r = ref;
		}
		else if (read == ',')
		{
			isforward = 0;
			r = ref;
		}
		else
			isforward = read == r ? 1 : 0;

		// find error allele
		errorbase=r;
		while (errorbase == r)
			errorbase = drawBase();

		if (errorbase == ref)
			read = isforward ? '.' : ',';
		else
			read = isforward ? errorbase : tolower(errorbase);
	}
}

int Assocsamp::minorid () const
{
	return _minorid;
}

double Assocsamp::minorfreq (int treatment)
{
	if (treatment <= _ntreatments)
		return _minorf[treatment];
	else
	{
		fprintf(stderr,"Invalid treatment %i passed to Assocsamp::minorfreq\n",treatment);
		_fail=1;
		return -1;
	}
}

void Assocsamp::printPileup (std::ostream& os, std::string chr, unsigned int pos, std::vector<Pileup>* indata)
{
	std::vector<Pileup>* pile = indata ? indata : &_pile;
	int i;

	if (pile->empty())
	{
		fprintf(stderr, "Empty vector of pileup information passed to Assocsamp::printPileup");
		return;
	}

	static std::vector<Pileup>::const_iterator iter;

	os << chr << "\t" << pos << "\t" << (*pile)[0].ref;

	for(iter=pile->begin(); iter!=pile->end(); ++iter)
	{
		if (iter->cov == 0)
			os << "\t" << "0" << "\t" << "*" << "\t" << "*";
		else
		{
			os << "\t" << iter->cov << "\t";

			for (i=0; i<iter->cov; ++i)
				os << iter->reads[i];

			os << "\t";

			for (i=0; i<iter->cov; ++i)
				os << iter->qual[i];
		}
	}

	os << "\n";
}

void Assocsamp::printParams (std::string* chr, unsigned int pos, std::ostream& os)
{
	os <<  *chr << "\t" << pos << "\t" << _fst << "\t" << _minorf[_ntreatments];
	for (int i=0; i<_ntreatments; ++i)
	{
		os << "\t" << _minorf[i];
	}
	os << "\n";
}
