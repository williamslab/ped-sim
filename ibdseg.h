// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
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
		    bool onlyGenetLen);
void printIBD(FILE *out, SimDetails &pedDetails, int fam,
	      vector< vector< vector<IBDRecord> > > *theSegs,
	      GeneticMap &map, bool sexSpecificMaps, bool onlyGenetLen);
void mergeSegments(vector<IBDRecord> &segs);
void printOneIBDSegment(FILE *out, SimDetails &pedDetails, int fam,
			int gen, int branch, int ind, IBDRecord &seg,
			int realStart, int realEnd, const char *type,
			GeneticMap &map, bool sexSpecificMaps,
			bool onlyGenetLen);
void clearTheSegs(SimDetails &pedDetails, 
		  vector< vector< vector<IBDRecord> > > *theSegs);

#endif // IBDSEG_H
