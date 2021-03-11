// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include <unordered_map>
#include "datastructs.h"
#include "geneticmap.h"
#include "fileorgz.h"

#ifndef BPVCFFAM_H
#define BPVCFFAM_H

using namespace std;

void readSexes(unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
	       uint32_t sexCount[2], const char *sexesFile);
template<typename O_TYPE = FILE *>
bool printSampleId(FILE *out, SimDetails &pedDetails, int fam, int gen,
		   int branch, int ind, bool printAllGens = false,
		   FileOrGZ<O_TYPE> *gzOut = NULL);
void printBPs(vector<SimDetails> &simDetails, Person *****theSamples,
	      GeneticMap &map, char *bpFile);
int printVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	     int totalFounderHaps, const char *inVCFfile, char *outFile,
	     GeneticMap &map, FILE *outs[2], vector<int> hapNumsBySex[2],
	     unordered_map<const char*,uint8_t,HashString,EqString> &sexes);
template<typename I_TYPE, typename O_TYPE>
int makeVCF(vector<SimDetails> &simDetails, Person *****theSamples,
	    int totalFounderHaps, const char *inVCFfile, char *outFileBuf,
	    GeneticMap &map, FILE *outs[2], vector<int> hapNumsBySex[2],
	    unordered_map<const char*,uint8_t,HashString,EqString> &sexes);
void getSampleIdsShuffHaps(vector<char*> &sampleIds,
		vector<uint8_t> &sampleSexes, vector<int> &shuffHaps,
		FILE *outs[2], vector<int> hapNumsBySex[2],
		unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
		char *&saveptr, int totalFounderHaps, const char *tab);
void printFam(vector<SimDetails> &simDetails, Person *****theSamples,
	      const char *famFile);
template<typename T>
void pop_front(vector<T> &vec);

#endif // BPVCFFAM_H
