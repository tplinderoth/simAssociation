/*
 * Argparse.h
 */

#ifndef ARGPARSE_H_
#define ARGPARSE_H_

#include <fstream>
#include <vector>

class Argparse
{
public:
	Argparse ();
	int parseInput (const int c, char** v, const char* version);
	void maininfo (const char* v);
	int fail();
	double freq; /* causative allele frequency */
	double minfreq; /* minimum allele frequency */
	double maxfreq; /* maximum allele frequency */
	double assoclevel; /* how many more times frequent the causative allele is in cases vs controls */
	std::vector<int> npools; /* number of pools for each treatment */
	std::vector<int> nind; /* diploid pool size for each treatment */
	double avgcov; /* average coverage per pool */
	double error; /* average sequencing error rate */
	double maxq; /* maximum quality score */
	double qoffset; /* quality score offset */
	unsigned int nsites; /* number of sites to simulate */
	std::string seqname; /* name of sequence, i.e. chromosome, scaffold, contig, etc. */
	int npops; /* number of populations to sample from */
	double fst; /* Fst when individuals from two populations are sampled*/
	std::string outfile; /* outfile name */
private:
	int _fail; /* flag for error status */
};

#endif /* ARGPARSE_H_ */
