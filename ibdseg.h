// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include <tuple>
#include "datastructs.h"
#include "geneticmap.h"

#ifndef IBDSEG_H
#define IBDSEG_H

using namespace std;

bool compInheritRecSamp(const InheritRecord &a, const InheritRecord &b);
bool compInheritRecStart(const InheritRecord &a, const InheritRecord &b);
bool compIBDRecord(const IBDRecord &a, const IBDRecord &b);
void locatePrintIBD(vector<SimDetails> &simDetails,
		    vector< vector< vector<InheritRecord> > > &hapCarriers,
		    GeneticMap &map, bool sexSpecificMaps, char *ibdFile,
		    vector< tuple<uint8_t,int,int,uint8_t,float> > *ibdSegs,
		    char *mrcaFile);
void printIBD(FILE *out, SimDetails &pedDetails, int fam,
	      vector< vector< vector<IBDRecord> > > *theSegs,
	      GeneticMap &map, bool sexSpecificMaps,
	      vector< tuple<uint8_t,int,int,uint8_t,float> > *ibdSegs,
	      FILE *mrcaOut);
void mergeSegments(vector<IBDRecord> &segs, bool retainFoundHap);
void printOneIBDSegment(FILE *out, SimDetails &pedDetails, int fam,
			int gen, int branch, int ind, IBDRecord &seg,
			int realStart, int realEnd, uint8_t ibdType,
			GeneticMap &map, bool sexSpecificMaps,
			vector<tuple<uint8_t,int,int,uint8_t,float> > *ibdSegs);
void printSegFounderId(FILE *mrcaOut, int foundHapNum, SimDetails &pedDetails,
		       int fam);
void clearTheSegs(SimDetails &pedDetails, 
		  vector< vector< vector<IBDRecord> > > *theSegs);

#endif // IBDSEG_H
