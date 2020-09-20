// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include <set>
#include "datastructs.h"

#ifndef READDEF_H
#define READDEF_H

using namespace std;

void readDef(vector<SimDetails> &simDetails, char *defFile);
void assignDefaultBranchParents(int prevGenNumBranches, int thisGenNumBranches,
				Parent **thisGenBranchParents, int prevGen,
				int *prevGenSpouseNum = NULL,
				vector<bool> *branchParentsAssigned = NULL);
bool readBranchSpec(int *numBranches, Parent **thisGenBranchParents,
		    int *thisGenNumSampsToPrint, int curGen,
		    int **sexConstraints, int **prevGenSpouseNum,
		    vector<bool> &branchParentsAssigned,
		    vector< set<Parent,ParentComp>* > &spouseDependencies,
		    const int i1Sex, const char *delim, char *&saveptr,
		    char *&endptr, int line);
void readParents(int *numBranches, int prevGen, int **sexConstraints,
		 int **prevGenSpouseNum,
		 vector< set<Parent,ParentComp>* > &spouseDependencies,
		 char *assignBranches, char *assignPar[2], Parent pars[2],
		 char *&fullAssignPar, const int i1Sex, char *&endptr,
		 int line);
void updateSexConstraints(int **sexConstraints, Parent pars[2],
			  int *numBranches,
			  vector< set<Parent,ParentComp>*> &spouseDependencies,
			  int line);
bool intersectNonEmpty(set<Parent,ParentComp> &a, set<Parent,ParentComp> &b);

#endif // READDEF_H
