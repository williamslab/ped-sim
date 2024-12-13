// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include "datastructs.h"
#include "geneticmap.h"
#include "fileorgz.h"

#ifndef BPVCFFAM_H
#define BPVCFFAM_H

using namespace std;

void readSexes(unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
               uint32_t sexCount[2], const char *sexesFile);
template<typename O_TYPE = FILE *>
bool printSampleId(FILE *out, SimDetails &pedDetails, int rep, int gen,
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
void assignSimId2Person(
    unordered_map<const char*,Person*,HashString,EqString> &simId2Person,
    vector<SimDetails> &simDetails, Person *****theSamples);
uint32_t setFoundersHapsToShuf(
          unordered_map<const char*,Person*,HashString,EqString> &simId2Person,
          unordered_map<const char*,uint32_t,HashString,EqString> &vcfId2index,
          unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
          vector<vector<int>> &vcfIdx2FounderHap,
          vector<int> hapNumsToShuf[3], vector<int> hapNumsBySex[2],
          int totalFounderHaps, uint32_t vcfSexCounts[3]);
uint32_t getSampleIdsShuffHaps(vector<char*> &vcfIds,
        vector<uint8_t> &vcfSampleSexes,
        unordered_map<const char*,Person*,HashString,EqString> &simId2Person,
        vector<vector<int>> &vcfIdx2FounderHap, FILE *outs[2],
        vector<int> hapNumsBySex[2],
        unordered_map<const char*,uint8_t,HashString,EqString> &sexes,
        char *&saveptr, int totalFounderHaps, const char *tab);
void printFam(vector<SimDetails> &simDetails, Person *****theSamples,
              const char *famFile);
template<typename T>
void pop_front(vector<T> &vec);

#endif // BPVCFFAM_H
