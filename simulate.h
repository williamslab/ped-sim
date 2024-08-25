// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include "datastructs.h"
#include "geneticmap.h"
#include "cointerfere.h"

#ifndef SIMULATE_H
#define SIMULATE_H

using namespace std;

// global random variables used in several other places
extern mt19937 randomGen;
extern uniform_int_distribution<int> coinFlip;
extern exponential_distribution<double> crossoverDist;

int simulate(vector<SimDetails> &simDetails, Person *****&theSamples,
	     GeneticMap &map, bool sexSpecificMaps, vector<COInterfere> &coIntf,
	     vector< vector< vector<InheritRecord> > > &hapCarriers,
	     vector<int> hapNumsBySex[2]);
void getPersonCounts(int curGen, int numGen, int branch, int **numSampsToPrint,
		     Parent **branchParents, int **branchNumSpouses,
		     int &numFounders, int &numNonFounders);
void generateHaplotype(Haplotype &toGenerate, Person &parent,
		       GeneticMap &map, vector<COInterfere> &coIntf,
		       unsigned int chrIdx,
		       vector< vector< vector<InheritRecord> > > &hapCarriers,
		       int ped, int rep, int curGen, int branch, int ind,
		       unsigned int fixedCOidxs[2]);
void copySegs(Haplotype &toGenerate, Person &parent, int &nextSegStart,
	      int switchPos, unsigned int curSegIdx[2], int &curHap,
	      unsigned int chrIdx,
	      vector< vector< vector<InheritRecord> > > &hapCarriers,
	      int ped, int rep, int curGen, int branch, int ind);
int getBranchNumSpouses(SimDetails &pedDetails, int gen, int branch);
void deleteTheSamples(vector<SimDetails> &simDetails, Person *****theSamples);


#endif // SIMULATE_H
