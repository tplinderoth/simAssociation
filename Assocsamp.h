/*
 * Assocsamp.h
 */

#ifndef ASSOCSAMP_H_
#define ASSOCSAMP_H_

#include <vector>
#include <string>
#include <fstream>

struct Pileup
{
	void allocate_read_storage (int depth);
	int cov; // coverage
	char ref; // reference allele
	char alt; // alterante allele
	std::string reads; // reads
	std::string qual; // quality scores
};

class Assocsamp
{
public:
	Assocsamp ();
	~Assocsamp ();
	int makePileup (double avgcov, double error, double maxqual, double offsetq, std::vector<Pileup>* pilevec=NULL); /* generate sequence data in pileup format from _npools */
	int quantitative_pool (int nind, int npools, double f [], int treatmentid, int ntreatments, int nsubpop, double fst=0.0);
	void genGenotypes (double f, int n, int geno []); /* expected [aa, aA, AA] genotype counts under HWE, f=p(a) */
	void treatGenos (double f, int n, double dominance, bool iscase, int geno []); /* expected [aa, aA, AA] counts for either a case or control sample */
	int& allele1_count (int pool, int treatmentid); /* return reference to count of allele1 of a specific pool and treatment */
	int& allele2_count (int pool, int treatmentid); /* return reference to count of allele2 of a specific pool and treatment */
	int allele1_count (int pool, int treatmentid) const; /* return count of allele1 of a specific pool and treatment */
	int allele2_count (int pool, int treatmentid) const; /* return count of allele2 of a specific pool and treatment */
	int minorid () const; /* return minor allele ID (allele 0 or allele 1) for a given treatment */
	double minorfreq (int treatment); /* return minor allele frequency for a given treatment */
	void sampleMultipops (int n, double fst, int npops, double f [], int totalgeno [], bool isancestral); /* sample genotypes from subpopulations */
	void setMinorfreq (); /* sets minor allele ID, _minorid, and frequencies, _minorf, based on the data stored in _pools */
	void printPileup (std::ostream& os, std::string chr, unsigned int pos, std::vector<Pileup>* indata=NULL); /* print pileup output */
	void printParams (std::string* chr, unsigned int pos, std::ostream& os); /* print simulation parameters */
private:
	void initialize_pool (int npool, int treatment);
	void initMinorfreq ();
	void initAlleleCounts ();
	void pileReads (double minorfreq, Pileup* pile); /* turn allele counts into string of pileup reads */
	void pileQuality (double rate, double maxq, double offset, Pileup* pile); /* generate phred scores for reads and create sequencing errors */
	void doError (int qscore, char ref, char alt, char& read); /* create sequencing error */
	char drawBase ();
	int arraysum(int a []);
	void setAlleleID (char* ref, char* alt);
	int _ntreatments; /* number of different treatments */
	int _npools; /* number of pools */
	std::vector<int**> _pools; /* each element/pool is [n causative alleles, n non-causative alleles]*/
	int** _counts; /* [total allele1 count, total allele2 count] for each treatment */
	int _minorid; /* minor allele (0 or 1) for site */
	double* _minorf; /* vector of minor allele frequency for treatment1, treatment2,..,treatmentn, site-wide */
	std::vector<Pileup> _pile; /* stores pileup version of the data in _pools */
	double _fst; /* fst between sampled populations */
	int _fail; /* flag indicating whether an error occurred */
};

#endif /* ASSOCSAMP_H_ */
