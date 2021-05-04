// ped-sim: pedigree simulation tool
//
// This program is distributed under the terms of the GNU General Public License

#include <vector>
#include <set>
#include <utility>
#include "datastructs.h"

#ifndef READDEF_H
#define READDEF_H

using namespace std;

void readDef(vector<SimDetails> &simDetails, char *defFile);
void finishLastDef(int numGen, SexConstraint **&sexConstraints,
		   vector< pair<set<Parent,ParentComp>,int8_t>* > &spouseDependencies);
void assignDefaultBranchParents(int prevGenNumBranches, int thisGenNumBranches,
				Parent **thisGenBranchParents, int prevGen,
				int *prevGenSpouseNum = NULL,
				vector<bool> *branchParentsAssigned = NULL);
bool readBranchSpec(int *numBranches, Parent **thisGenBranchParents,
		    int *thisGenNumSampsToPrint, int curGen,
		    SexConstraint **sexConstraints, int **prevGenSpouseNum,
		    vector<bool> &branchParentsAssigned,
		    vector< pair<set<Parent,ParentComp>,int8_t>* > &spouseDependencies,
		    const int i1Sex, const char *delim, char *&saveptr,
		    char *&endptr, int line);
void assignBranch(bool parentAssign, bool noPrint, int sexToAssign, int curGen,
		  int branch, SexConstraint **sexConstraints,
		  Parent **thisGenBranchParents, int *thisGenNumSampsToPrint,
		  int thisGenNumBranches, vector<bool> &branchParentsAssigned,
		  Parent pars[2], int line, bool &warningGiven);
void readParents(int *numBranches, int prevGen, SexConstraint **sexConstraints,
		 int **prevGenSpouseNum,
		 vector< pair<set<Parent,ParentComp>,int8_t>* > &spouseDependencies,
		 char *assignBranches, char *assignPar[2], Parent pars[2],
		 char *&fullAssignPar, const int i1Sex, char *&endptr,
		 int line);
void updateSexConstraints(SexConstraint **sexConstraints, Parent pars[2],
			  int *numBranches,
			  vector< pair<set<Parent,ParentComp>,int8_t>* > &spouseDependencies,
			  int line);
bool intersectNonEmpty(set<Parent,ParentComp> &a, set<Parent,ParentComp> &b);
void initSexConstraints(SexConstraint *newSexConstraints, int numBranches);

#endif // READDEF_H
