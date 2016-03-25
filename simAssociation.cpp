/*
 * simAssociation.cpp
 */

#include "simAssociation.h"
#include "Assocsamp.h"
#include "generalUtils.h"
#include "randomVar.h"

int main (int argc, char** argv)
{
	int status=0;
	Argparse params;

	status=params.parseInput(argc,argv,version);
	if (status)
		return status;

	status=simPileup(&params);
	if (status)
	{
		fprintf(stderr,"Program terminated with errors\n");
		return status;
	}

	fprintf(stderr,"finished!\n");

	return status;
}

int simPileup (Argparse* par)
{
	Assocsamp sample;
	std::fstream pile_fs;
	std::ostream pile_os(std::cout.rdbuf());
	std::fstream param_fs;
	std::ostream param_os(std::cout.rdbuf());
	int ntreat = par->nind.size();
	int numpops=par->npops;
	double fst=par->fst;
	unsigned int i;
	int j,k;
	int rc = 0;
	double f [ntreat][numpops];

	// set output streams

	std::string paramout = par->outfile;
	paramout += ".param";
	if (outStreams(par->outfile, pile_os, pile_fs))
		return -1;
	if(outStreams(paramout, param_os, param_fs))
		return -1;

	// simulate pileup data

	for (i=0; i<par->nsites; ++i)
	{
		// set allele frequencies
		for(j=0; j<ntreat; ++j)
		{
			for (k=0; k<numpops; ++k)
			{
				if (j == 0)
				{
					f[j][k] = par->freq == 82.0 ? decimalUnifBound(par->minfreq, par->maxfreq) : par->freq;
					if (numpops > 1 && par->fst != 0.0)
						f[j][k] = randomVar::rbeta((((1.0-fst)/fst)*f[j][k]), (((1.0-fst)/fst)*(1.0-f[j][k])));
				}
				else
					f[j][k] = (f[j][k]=f[0][k]*par->assoclevel) > 1.0 ? 1.0 : f[j][k];
			}
		}

		// sample individuals
		for (j=0; j<ntreat; ++j)
		{
			if(sample.quantitative_pool(par->nind[j], par->npools[j], f[j], j, ntreat, par->npops, fst))
			{
				fprintf(stderr,"Unable to carry out Assocsamp::quantitative_pool\n");
				rc = -1;
				break;
			}
		}

		// generate pileup data from sampled genotypes
		if (sample.makePileup(par->avgcov, par->error, par->maxq, par->qoffset))
		{
			fprintf(stderr,"Unable to make pileup data\n");
			rc = -1;
			break;
		}
		sample.printPileup(pile_os, par->seqname, i+1);
		sample.printParams(&par->seqname, i+1, param_os);
	}

	// close files
	pile_fs.close();
	param_fs.close();

	return rc;
}

int outStreams (const std::string outfile, std::ostream& os, std::fstream& fs)
{
	// open output steams
	if (!outfile.empty())
	{
		if (getFILE(fs, outfile.c_str(), "out"))
		{
			os.rdbuf(fs.rdbuf()); // switch output stream's buffer to output file's buffer
			std::cerr << "Dumping results to " << outfile << "\n";
		}
		else
		{
			fprintf(stderr,"Problem opening output file ...\n");
			return -1;
		}
	}
	else
	{
		fprintf(stderr,"No output file name passed to outStreams\n");
		return -1;
	}

	return 0;
}
