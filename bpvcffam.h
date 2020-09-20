// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include "datastructs.h"
#include "geneticmap.h"
#include "fileorgz.h"

#ifndef BPVCFFAM_H
#define BPVCFFAM_H

using namespace std;

template<typename IO_TYPE = FILE *>
bool printSampleId(FILE *out, SimDetails &pedDetails, int fam, int gen,
		   int branch, int ind, bool printAllGens = false,
		   FileOrGZ<IO_TYPE> *gzOut = NULL);
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      GeneticMap &map, char *bpFile);
void printVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	      int totalFounderHaps, const char *inVCFfile, char *outFile,
	      GeneticMap &map, FILE *outs[2]);
template<typename IO_TYPE>
void makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, const char *inVCFfile, char *outFileBuf,
	     GeneticMap &map, FILE *outs[2]);
void printFam(vector<SimDetails> &simDetails, Person *****theSamples,
	      char *famFile);
template<typename T>
void pop_front(vector<T> &vec);

#endif // BPVCFFAM_H
