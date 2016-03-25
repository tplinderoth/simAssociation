/*
 * Argparse.cpp
 */

#include "Argparse.h"
#include "generalUtils.h"
#include <sstream>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <math.h>

Argparse::Argparse ()
	: freq (82.0),
	  minfreq(0.02),
	  maxfreq(0.98),
	  assoclevel(1.0),
	  avgcov(5.0),
	  error(0.001),
	  maxq(41.0),
	  qoffset(33.0),
	  nsites(1),
	  seqname("chr1"),
	  npops(1),
	  fst(0.0),
	  _fail(0)
{}

int Argparse::parseInput (const int c, char** v, const char* version)
{
	int argPos = 1;
	int increment = 0;
	int i;

	if (c < 2)
	{
		maininfo(version);
		return 1;
	}

	while (argPos < c)
	{
		if (strcmp(v[argPos],"-help") == 0)
		{
			maininfo(version);
			return 1;
		}
		else if (strcmp(v[argPos], "-outfile") == 0)
		{
			outfile = v[argPos+1];
		}
		else if (strcmp(v[argPos], "-npools") == 0)
		{
			npools.reserve(2);
			i=1;
			while (argPos+i < c && isdigit(v[argPos+i][0]))
			{
				if (atoi(v[argPos+i]) < 1)
				{
					fprintf(stderr,"Treatment %i has no samples\n",i);
					return -1;
				}
				npools.push_back(atoi(v[argPos+i]));
				++i;
			}
			increment=i-2;
		}
		else if (strcmp(v[argPos], "-nind") == 0)
		{
			nind.reserve(2);
			i=1;
			while (argPos+i < c && isdigit(v[argPos+i][0]))
			{
				if (atoi(v[argPos+i]) < 1)
				{
					fprintf(stderr,"Treatment %i has pool size of zero\n",i);
					return -1;
				}
				nind.push_back(atoi(v[argPos+i]));
				++i;
			}
			increment=i-2;
		}
		else if (strcmp(v[argPos], "-avgcov") == 0)
		{
			avgcov = atof(v[argPos+1]);
			if (avgcov <= 0.0)
			{
				fprintf(stderr, "average coverage must be > 0.0\n");
				return -1;
			}

		}
		else if (strcmp(v[argPos], "-error")==0)
		{
			error=atof(v[argPos+1]);
			if (error < 0.0 || error > 1.0)
			{
				fprintf(stderr,"Sequencing error rate must be in [0,1]\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-freq") == 0)
		{
			if (strcmp(v[argPos+1], "R") == 0 || strcmp(v[argPos+1], "r") == 0 )
			{
				freq = 82.0; // R ascii decimal value
			}
			else
			{
				freq = atof(v[argPos+1]);
				if (freq < 0.0 || freq > 1.0)
				{
					fprintf(stderr,"Invalid allele frequency %f\n",freq);
					return -1;
				}
			}
		}
		else if (strcmp(v[argPos], "-minfreq") == 0)
		{
			minfreq=atof(v[argPos+1]);

		}
		else if (strcmp(v[argPos],"-maxfreq")==0)
		{
			maxfreq=atof(v[argPos+1]);
		}
		else if (strcmp(v[argPos], "-assoclevel") == 0)
		{
			assoclevel=atof(v[argPos+1]);
			if (assoclevel < 0.0 )
			{
				fprintf(stderr, "-assoclevel argument must be >= 0\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-maxq") == 0)
		{
			maxq=atof(v[argPos+1]);
			if (maxq < 0.0 )
			{
				fprintf(stderr,"-maxq argument must be >= 0.0\n");
				return -1;
			}

		}
		else if (strcmp(v[argPos], "-qoffset") == 0)
		{
			qoffset=atof(v[argPos+1]);
		}
		else if (strcmp(v[argPos], "-nsites") == 0)
		{
			std::stringstream strvalue;
			strvalue << v[argPos+1];
			strvalue >> nsites;
			if (nsites < 0)
			{
				fprintf(stderr,"Invalid number of sites to simulate: %u",nsites);
				return -1;
			}
		}
		else if (strcmp(v[argPos], "-seqname") == 0)
		{
			seqname = v[argPos+1];
			if (seqname.empty())
			{
				fprintf(stderr,"No sequence name supplied\n");
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-fst")==0)
		{
			fst = atof(v[argPos+1]);
			if (fst < 0.0 || fst > 1.0)
			{
				fprintf(stderr,"Invalid Fst value %f",fst);
				return -1;
			}
		}
		else if (strcmp(v[argPos],"-npops")==0)
		{
			npops = atoi(v[argPos+1]);
			if (npops < 1)
			{
				fprintf(stderr,"number of populations to sample must be at least 1\n");
				return -1;
			}
		}
		else
		{
			fprintf(stderr, "Unknown option: %s\n", v[argPos]);
			return -1;
		}
		argPos += 2 + increment;
		increment = 0;
	}

	if (-10*log10(error) > maxq)
	{
		fprintf(stderr, "Phred-scaled error rate %f exceeds max base quality %f",error,maxq);
		return -1;
	}

	if (qoffset < 0.0 || qoffset >= maxq)
	{
		fprintf(stderr,"Invalid quality score offset _qoffset %f",qoffset);
		return -1;
	}

	if (minfreq < 0.0 || maxfreq > 1.0 || maxfreq < minfreq)
	{
		fprintf(stderr,"Invalid bounds on causative allele frequency [%f %f]\n",minfreq,maxfreq);
		return-1;
	}

	return 0;
}

int Argparse::fail()
{
	return _fail;
}

void Argparse::maininfo (const char* v)
{
	int w = 12;
	std::cerr << "\nsimAssociation version " << v << "\n\nUsage: simAssociation [options]\n"
	<< "\nInput:\n"
	<< "\n" << std::setw(w) << std::left << "-outfile" << std::setw(w) << "FILE" << "File to output results to"
	<< "\n" << std::setw(w) << std::left << "-freq" << std::setw(w) << "FLOAT|R" << "Frequency of causative allele; 'R' specifies random frequency bounded by -minfreq and -maxfreq [" << static_cast<char>(freq) << "]"
	<< "\n" << std::setw(w) << std::left << "-minfreq" << std::setw(w) << "FLOAT" << "Minimum causative allele frequency [" << minfreq << "]"
	<< "\n" << std::setw(w) << std::left << "-maxfreq" << std::setw(w) << "FLOAT" << "Maximum causative allele frequency [" << maxfreq << "]"
	<< "\n" << std::setw(w) << std::left << "-assoclevel" << std::setw(w) << "FLOAT" << "How many more times frequent the causative allele is in the case population [" << assoclevel << "]"
	<< "\n" << std::setw(w) << std::left << "-npops" << std::setw(w) << "INT" << "Sample an equal number of individuals from among INT populations [" << npops << "]"
	<< "\n" << std::setw(w) << std::left << "-npools" << std::setw(w) << "INT" << "An array of the number of pools for treatment1, treatment2,...,treatmentn"
	<< "\n" << std::setw(w) << std::left << "-nind" << std::setw(w) << "INT" << "An array of the diploid pools size for treatment1, treatment2,...,treatmentn"
	<< "\n" << std::setw(w) << std::left << "-avgcov" << std::setw(w) << "FLOAT" << "Average sequencing depth per pool [" << avgcov << "]"
	<< "\n" << std::setw(w) << std::left << "-error" << std::setw(w) << "FLOAT" << "Average sequencing error rate [" << error << "]"
	<< "\n" << std::setw(w) << std::left << "-maxq" << std::setw(w) << "FLOAT" << "Maximum Phred-scaled quality score [" << maxq << "]"
	<< "\n" << std::setw(w) << std::left << "-qoffset" << std::setw(w) << "FLOAT" << "Minimum possible ASCII decimal value used to encode quality scores [" << qoffset << "]"
	<< "\n" << std::setw(w) << std::left << "-nsites" << std::setw(w) << "INT" << "Number of sites to simulate [" << nsites << "]"
	<< "\n" << std::setw(w) << std::left << "-seqname" << std::setw(w) << "FLOAT" << "Name of sequenced chromosome, scaffold, contig, etc. [" << seqname << "]"
	<< "\n" << std::setw(w) << std::left << "-fst" << std::setw(w) << "FLOAT" << "Fst between sampled populations [" << fst << "]"
	<< "\n\nOutput parameter file format:"
	<< "\n(1) chromosome name"
	<< "\n(2) site number"
	<< "\n(3) Fst"
	<< "\n(4) sample-wide minor allele frequency"
	<< "\n(5+) treatment-specific minor allele frequency"
	<< "\n\n";
}
