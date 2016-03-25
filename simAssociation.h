/*
 * simAssociation.h
 */

#ifndef SIMASSOCIATION_H_
#define SIMASSOCIATION_H_

#include "Argparse.h"

const char* version = "0.0.3";

int simPileup (Argparse* par);
int outStreams (const std::string outfile, std::ostream& os, std::fstream& fs);

#endif /* SIMASSOCIATION_H_ */
